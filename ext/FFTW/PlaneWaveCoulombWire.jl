struct PlaneWaveCoulombWire{dctT<:FFTW.AbstractFFTs.Plan{Float64}} <: PlaneWaveCoulomb
    lattice::Lattice{Float64}
    period::Vec3{Bool}
    Glist::Vector{ReducedCoordinates{Int}}
    fft_grid::NTuple{2,Int}
    α::Float64
    zindex::Int
    xyindex::SVector{2,Int}
    Rx::Float64
    Ry::Float64
    R::Float64
    rlatticez::Float64
    Gz::Vector{Int}
    Gxy::Vector{SVector{2,Int}}
    Gzcar::Vector{Float64}
    dct_K0::dctT
    r::Matrix{Float64}
    coff::Float64
end
function (pwc::PlaneWaveCoulombWire)(q::AbstractVector{<:Real}, Gidx::AbstractVector{<:Integer})
    coulomb = Vector{Float64}(undef, length(Gidx))
    return pwc(coulomb, q, Gidx)
end
@inbounds function (pwc::PlaneWaveCoulombWire)(coulomb::AbstractVector{<:Real}, q::AbstractVector{<:Real}, Gidx::AbstractVector{<:Integer})
    qcar = pwc.rlatticez * q[pwc.zindex]
    Gz_q = view(pwc.Gz, Gidx)
    Gxy_q = view(pwc.Gxy, Gidx)
    Gzcar_q = view(pwc.Gzcar, Gidx)

    # unique Gz, and its index in Glist_q
    Gz_unique, Gz_indices = unique_groups(Gz_q)

    Threads.@threads for iGz in eachindex(Gz_unique)
        realgrid = Matrix{Float64}(undef, pwc.fft_grid[1], pwc.fft_grid[2])
        iGs_iGz = Gz_indices[iGz]
        qGz = qcar + Gzcar_q[iGs_iGz[1]]
        realgrid .= besselk.(0, qGz .* pwc.r)
        pwc.dct_K0 * realgrid
        for iG in iGs_iGz
            Gxy = Gxy_q[iG]
            #TODO: 这里的G̃x有可能造成越界访问，需要周期映射回来
            G̃x, G̃y = abs(Gxy[1]) + 1, abs(Gxy[2]) + 1
            coulomb[iG] = realgrid[G̃x, G̃y] * pwc.coff
        end
    end

    # iGz_chunks = _split_tasks(length(Gz_unique), pwc.nchunk)
    # Threads.@threads for ichunk in Base.OneTo(pwc.nchunk)
    # 	realgrid = pwc.realgrids[ichunk]
    # 	for iGz in iGz_chunks[ichunk]
    # 		iGs_iGz = Gz_indices[iGz]
    # 		qGz = qcar + Gzcar_q[iGs_iGz[1]]
    # 		realgrid .= besselk.(0, qGz .* pwc.r)
    # 		pwc.dct_K0 * realgrid
    # 		qGz2 = qGz * qGz
    # 		for iG in iGs_iGz
    # 			Gxy = Gxy_q[iG]
    # 			G̃x, G̃y = abs(Gxy[1]) + 1, abs(Gxy[2]) + 1
    # 			coulomb[iG] = realgrid[G̃x, G̃y] * pwc.coff
    # 		end
    # 	end
    # end

    return coulomb
end
function (pwc::PlaneWaveCoulombWire)(::Val{:head}, qgrid::Union{AbstractVector,Tuple})
    R² = pwc.R * pwc.R
    nq = qgrid[pwc.zindex]
    qsz = pwc.rlatticez / nq / 2

    function integrand(xyq)
        r² = xyq[1] * xyq[1] + xyq[2] * xyq[2]
        r² < R² && return 0.0
        return besselk(0, xyq[3] * sqrt(r²))
    end
    (result, err) = hcubature(integrand, (0, 0, 0), (pwc.Rx/2, pwc.Ry/2, qsz); norm=abs, rtol=1e-8, atol=1e-10)

    function integrand_cylinder(k)
        kR = k * pwc.R
        return (1 - kR * besselk(1, kR)) / (k * k)
    end
    (result_cylinder, err) = quadgk(integrand_cylinder, 0, qsz; rtol=1e-8, atol=1e-10)

    average = (result * 2 / π + result_cylinder) / qsz
    return pwc.coff * average
end
function (pwc::PlaneWaveCoulombWire)(::Val{:head}, qgrid::Union{AbstractVector,Tuple}, ϵ_∞::Real, normq0::Real)
    R² = pwc.R * pwc.R
    nq = qgrid[pwc.zindex]
    qsz = pwc.rlatticez / nq / 2
    function integrand_vq(q, r)
        r² = r[1] * r[1] + r[2] * r[2]
        r² < R² && return 0.0
        return besselk(0, q * sqrt(r²)) # requier positive q0
    end
    half_Rx = pwc.Rx / 2
    half_Ry = pwc.Ry / 2
    function vq(q)
        (result_q, err) = hcubature(r->integrand_vq(q, r), (0, 0), (half_Rx, half_Ry); norm=abs, rtol=1e-11, atol=1e-12)
        qR = q * pwc.R
        result_circle_q = (1 - qR * besselk(1, qR)) / (q * q)
        return result_q * 2 / π + result_circle_q
    end
    v_q0 = vq(normq0)
    factor_q0 = v_q0 * normq0 * normq0
    γ = (1 / ϵ_∞ - 1) * exp(pwc.α * normq0) / factor_q0

    function integrand(q)
        vq_value = vq(q)
        ϵq_value = 1 + γ * exp(-head.α * q) * q * q * vq_value
        return vq_value / ϵq_value
    end

    (result, err) = quadgk(integrand, 0, qsz; norm=abs, rtol=1e-8, atol=1e-9)
    average = result / qsz
    return pwc.coff * average
end
"""
    PlaneWaveCoulomb(:wire, lattice::Union{Lattice,AbstractMatrix{<:Real}}, Glist; period=[false, false, true], fft_grid, α::Real=0)

Construct a PlaneWaveCoulomb object for 1D systems with wire truncation.
"""
function PlaneWaveCoulomb(::Val{:wire}, lattice::Union{Lattice,AbstractMatrix{<:Real}}, Glist;
    period=[false, false, true], fft_grid, α::Real=0)

    check_lattice_period(lattice, period)
    period = Vec3{Bool}(period)
    @assert count(period) == 1 "Wrong period, wire truncation is used in 1D system."
    if period[1]
        zindex = 1
        xyindex = SVector{2,Int}(2, 3)
    elseif period[2]
        zindex = 2
        xyindex = SVector{2,Int}(3, 1)
    else
        zindex = 3
        xyindex = SVector{2,Int}(1, 2)
    end

    lattice = Lattice(parent(lattice))
    rlattice = reciprocal(lattice)
    Rx = lattice[xyindex[1], xyindex[1]]
    Ry = lattice[xyindex[2], xyindex[2]]
    R = min(Rx, Ry) / 2

    Ω = abs(det(parent(lattice)))

    rlatticez = rlattice[zindex, zindex]

    Gz = map(Glist) do G
        G[zindex]
    end
    Gxy = map(Glist) do G
        G[xyindex]
    end
    Gzcar = Gz .* rlatticez

    half_lx = lattice[xyindex[1], xyindex[1]] / 2
    half_ly = lattice[xyindex[2], xyindex[2]] / 2

    fft_grid = (Int(fft_grid[xyindex[1]]), Int(fft_grid[xyindex[2]]))
    x = [half_lx * (i-0.5) / fft_grid[1] for i in 1:fft_grid[1]]
    y = [half_ly * (i-0.5) / fft_grid[2] for i in 1:fft_grid[2]]
    r = [sqrt(xx^2 + yy^2) for xx in x, yy in y]

    coff = half_lx * half_ly / (fft_grid[1] * fft_grid[2]) * 2 * CoulombScale / Ω

    FFTW_nt = FFTW.get_num_threads()
    FFTW.set_num_threads(1)
    # plan会保持建立时的线程数运行。
    dct_K0 = FFTW.plan_r2r!(Matrix{Float64}(undef, fft_grid[1], fft_grid[2]), FFTW.REDFT10)
    FFTW.set_num_threads(FFTW_nt)

    return PlaneWaveCoulombWire(lattice, period, Glist, fft_grid, α,
        zindex, xyindex, Rx, Ry, R, rlatticez, Gz, Gxy, Gzcar, dct_K0, r, coff)
end
