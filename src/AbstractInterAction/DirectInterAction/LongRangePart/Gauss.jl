function LongRangePart(::Val{:gauss}, lattice::Lattice,
    orblocat::Vector{<:ReducedCoordinates}; norb::Integer=length(orblocat),
    α::Real=1, αmat::AbstractMatrix{<:Real}=fill(α, (norb, norb)),
    ϵ::Real=1, ϵmat::AbstractMatrix{<:Real}=fill(ϵ, (norb, norb)),
    δ::Real=1e-8)

    Ω = abs(det(parent(lattice)))

    rlattice = reciprocal(lattice)
    rldot = parent(rlattice)
    rldot = Mat3(transpose(rldot) * rldot)

    αmat = _LongRangePartGauss_symmetrize_αmat(αmat)
    φ = [ReciprocalCoulomb(Val(:gauss); ϵ=ϵmat[i, j], α=αmat[i, j], Ω) for i in 1:norb, j in 1:norb]

    Gnorm_max = sqrt(-4 * maximum(αmat)^2 * log(δ)) + norm(parent(rlattice))
    G_frac = gridindex(rlattice, Gnorm_max, [true, true, true])
    nG = length(G_frac)

    G0idx = findfirst(iszero, G_frac)
    Gn0idx = setdiff(1:nG, G0idx)
    nGminus1 = nG - 1
    G_frac = G_frac[Gn0idx]
    G_car = map(G -> rlattice * G, G_frac)
    G_norm2 = map(G -> G[1] * G[1] + G[2] * G[2] + G[3] * G[3], G_car)

    Δorb = Matrix{ReducedCoordinates{Float64}}(undef, norb, norb)
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

    return LongRangePart3D{eltype(φ)}(lattice, rlattice, rldot, Ω, norb, orblocat, φ,
        nG, nGminus1, G_frac, G_car, G_norm2, Δorb, G_orb_phase)
end
function _LongRangePartGauss_symmetrize_αmat(αmat)
    norb = size(αmat, 1)
    αmat_sym = Matrix{Float64}(undef, norb, norb)
    for j in 2:norb, i in 1:(j-1)
        αmat_sym[i, j] = αmat[i, j]
        αmat_sym[j, i] = αmat[i, j]
    end
    for i in 1:norb
        αmat_sym[i, i] = αmat[i, i]
    end
    return αmat_sym
end
