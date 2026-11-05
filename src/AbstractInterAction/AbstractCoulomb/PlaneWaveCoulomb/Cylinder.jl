struct PlaneWaveCoulombCylinder <: PlaneWaveCoulomb
    lattice::Lattice{Float64}
    period::Vec3{Bool}
    Glist::Vector{ReducedCoordinates{Int}}
    R::Float64
    α::Float64
    zindex::Int
    xyindex::SVector{2,Int}
    rlatticez::Float64
    Gz::Vector{Int}
    Gxy2::Vector{Float64}
    GRJ1GR::Vector{Float64}
    J0GR::Vector{Float64}
    Ω::Float64
    coff::Float64
end
function (pwc::PlaneWaveCoulombCylinder)(q::AbstractVector{<:Real}, Gidx::AbstractVector{<:Integer})
    coulomb = Vector{Float64}(undef, length(Gidx))
    return pwc(coulomb, q, Gidx)
end
function (pwc::PlaneWaveCoulombCylinder)(coulomb::AbstractVector{<:Real}, q::AbstractVector{<:Real}, Gidx::AbstractVector{<:Integer})
    qz = q[pwc.zindex]
    Gz_q = view(pwc.Gz, Gidx)
    Gxy2_q = view(pwc.Gxy2, Gidx)
    GRJ1GR_q = view(pwc.GRJ1GR, Gidx)
    J0GR_q = view(pwc.J0GR, Gidx)
    Threads.@threads for iG in eachindex(Gidx)
        qG_z = pwc.rlatticez * (qz + Gz_q[iG])
        qG2 = qG_z * qG_z + Gxy2_q[iG]
        qG_z_R = abs(qG_z) * pwc.R
        factor = 1.0 + GRJ1GR_q[iG] * besselk(0, qG_z_R) - qG_z_R * J0GR_q[iG] * besselk(1, qG_z_R)
        coulomb[iG] = pwc.coff * factor / qG2
    end
    return coulomb
end
function (pwc::PlaneWaveCoulombCylinder)(::Val{:head}, qgrid::Union{AbstractVector,Tuple})

    nq = qgrid[pwc.zindex]
    qsz = pwc.rlatticez / nq / 2

	# TODO: 这里的积分有可能产生Inf，如果其采样点采集到了十分接近0的位置
    function integrand(k)
        kR = k * pwc.R
        return (1 - kR * besselk(1, kR)) / (k * k)
    end
    (result, err) = quadgk(integrand, 0, qsz; rtol=1e-8, atol=1e-10)

    # 被积函数为偶函数(省略了取模)，积分区域只积一半[0, qsz]，积分区域的长度也只需要qsz
    average = result / qsz
    return pwc.coff * average
end
function (pwc::PlaneWaveCoulombCylinder)(::Val{:head}, qgrid::Union{AbstractVector,Tuple}, ϵ_∞::Real, normq0::Real)
    error("Don't support cylinder head with screening.")
    # TODO
end
"""
    PlaneWaveCoulomb(:cylinder, lattice::Union{Lattice,AbstractMatrix{<:Real}}, Glist; period=[false, false, true], R::Union{Real,Nothing}=nothing, α::Real=0)

Construct a PlaneWaveCoulomb object for 1D systems with cylinder truncation.
"""
function PlaneWaveCoulomb(::Val{:cylinder}, lattice::Union{Lattice,AbstractMatrix{<:Real}}, Glist;
    period=[false, false, true], R::Union{Real,Nothing}=nothing, α::Real=0)

    check_lattice_period(lattice, period)
    period = Vec3{Bool}(period)
    @assert count(period) == 1 "Wrong period, cylinder truncation is used in 1D system."
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

    rlatticez = rlattice[zindex, zindex]
    rlattice_perp = SMatrix{2,2,Float64,4}(rlattice[xyindex, xyindex])
    rldot_perp = transpose(rlattice_perp) * rlattice_perp

    if isnothing(R)
        lattice_perp = lattice[xyindex, xyindex]
        h₁ = lattice_perp[1, 1]
        h₂ = lattice_perp[2, 2]
        R = min(h₁, h₂) / 2
    end

    Glist = ReducedCoordinates{Int}.(Glist)
    Gz = map(Glist) do G
        G[zindex]
    end

    nGlist = length(Glist)
    Gxy2 = Vector{Float64}(undef, nGlist)
    GRJ1GR = Vector{Float64}(undef, nGlist)
    J0GR = Vector{Float64}(undef, nGlist)
    Threads.@threads for iG in Base.OneTo(nGlist)
        Gxy = Glist[iG][xyindex]
        Gxy2[iG] = Gxy ⋅ (rldot_perp * Gxy)
        GR = sqrt(Gxy2[iG]) * R
        GRJ1GR[iG] = GR * besselj1(GR)
        J0GR[iG] = besselj0(GR)
    end

    Ω = abs(det(parent(lattice)))
    coff = CoulombScale * 4π / Ω

    return PlaneWaveCoulombCylinder(lattice, period, Glist, R, α,
        zindex, xyindex, rlatticez, Gz, Gxy2, GRJ1GR, J0GR, Ω, coff)
end
