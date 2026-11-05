function LongRangePart(::Val{:keldysh_gauss}, lattice::Lattice,
    orblocat::Vector{<:ReducedCoordinates}; norb::Integer=length(orblocat),
    period=Bool[1, 1, 0],
    α::Real=1, αmat::AbstractMatrix{<:Real}=fill(α, (norb, norb)),
    αz::Real=1, αzmat::AbstractMatrix{<:Real}=fill(αz, (norb, norb)),
    r₀::Real=0, r₀mat::AbstractMatrix{<:Real}=fill(r₀, (norb, norb)),
    ϵ₁::Real=1, ϵ₂::Real=1, δ::Real=1e-8)

    period = Vec3{Bool}(period)
    @assert count(period) == 2 "Exactly two dimensions should be periodic."
    if !period[1]
        xyindex = SVector{2,Int}(2, 3)
        zindex = 1
    elseif !period[2]
        xyindex = SVector{2,Int}(3, 1)
        zindex = 2
    else
        xyindex = SVector{2,Int}(1, 2)
        zindex = 3
    end

    @assert all(iszero, lattice[xyindex, zindex]) "Wrong lattice: the non-periodic direction should be perpendicular to the periodic plane."
    @assert all(iszero, lattice[zindex, xyindex]) "Wrong lattice: the non-periodic direction should be perpendicular to the periodic plane."

    rlattice = reciprocal(lattice)
    lattice2D = SMatrix{2,2,Float64,4}(lattice[xyindex, xyindex])
    rlattice2D = SMatrix{2,2,Float64,4}(rlattice[xyindex, xyindex])
    rldot = SMatrix{2,2,Float64,4}(transpose(rlattice2D) * rlattice2D)

    𝐚 = lattice2D[:, 1]
    𝐛 = lattice2D[:, 2]
    S = abs(𝐚[1] * 𝐛[2] - 𝐚[2] * 𝐛[1])

    orblocat = map(x->SVector{2,Float64}(x[xyindex]), orblocat)

    αmat, αzmat = _LongRangePartKeldyshGauss_symmetrize_αmat(αmat, αzmat)
    φ = [ReciprocalCoulomb(:keldysh_gauss; r₀=r₀mat[i, j], ϵ₁, ϵ₂, α=αmat[i, j], αz=αzmat[i, j], S) for i in 1:norb, j in 1:norb]

    Gnorm_max = sqrt(-4 * maximum(αmat)^2 * log(δ)) + norm(rlattice2D)
    period = [false, false, false]
    period[xyindex] .= true
    G_frac = gridindex(rlattice, Gnorm_max, period)
    G_frac = map(G -> SVector{2,Int}(G[xyindex]), G_frac)
    nG = length(G_frac)

    G0idx = findfirst(iszero, G_frac)
    Gn0idx = setdiff(1:nG, G0idx)
    nGminus1 = nG - 1
    G_frac = G_frac[Gn0idx]
    G_car = map(G -> rlattice2D * G, G_frac)
    G_norm2 = map(G -> G[1] * G[1] + G[2] * G[2], G_car)

    Δorb = Matrix{SVector{2,Float64}}(undef, norb, norb)
    G_orb_phase = Array{ComplexF64}(undef, nGminus1, norb, norb)
    for i in 1:norb, j in 1:i
        dorb = orblocat[i] - orblocat[j]
        Δorb[i, j] = dorb
        Δorb[j, i] = -dorb
        for (Gi, G) in enumerate(G_frac)
            G_orb_phase[Gi, i, j] = cispi(2 * (G ⋅ Δorb[i, j]))
        end
        G_orb_phase[:, j, i] .= conj.(G_orb_phase[:, i, j])
    end

    return LongRangePart2D{eltype(φ)}(lattice, rlattice, period, xyindex, zindex,
        lattice2D, rlattice2D, rldot, S, norb, orblocat, φ,
        nG, nGminus1, G_frac, G_car, G_norm2, Δorb, G_orb_phase)
end
function _LongRangePartKeldyshGauss_symmetrize_αmat(αmat, αzmat)
    norb = size(αmat, 1)
    αmat_sym = Matrix{Float64}(undef, norb, norb)
    αzmat_sym = Matrix{Float64}(undef, norb, norb)
    for j in 2:norb, i in 1:(j-1)
        αmat_sym[i, j] = αmat[i, j]
        αmat_sym[j, i] = αmat[i, j]
        αzmat_sym[i, j] = αzmat[i, j]
        αzmat_sym[j, i] = αzmat[i, j]
    end
    for i in 1:norb
        αmat_sym[i, i] = αmat[i, i]
        αzmat_sym[i, i] = αzmat[i, i]
    end
    return αmat_sym, αzmat_sym
end
