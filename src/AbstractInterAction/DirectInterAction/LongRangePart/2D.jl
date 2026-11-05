struct LongRangePart2D{φT<:ReciprocalCoulomb} <: LongRangePart
    # wannier structure.
    lattice::Lattice{Float64}
    rlattice::ReciprocalLattice{Float64}
    period::Vec3{Bool}
    xyindex::SVector{2,Int}
    zindex::Int
    lattice2D::SMatrix{2,2,Float64,4}
    rlattice2D::SMatrix{2,2,Float64,4}
    rldot::SMatrix{2,2,Float64,4}
    S::Float64
    norb::Int
    orblocat::Vector{SVector{2,Float64}}
    φ::Matrix{φT}
    # G grid for the summation of long-range correction.
    nG::Int
    # only contain nonzero G vectors, the length equal to nG-1.
    nGminus1::Int
    G_frac::Vector{SVector{2,Int}}
    G_car::Vector{SVector{2,Float64}}
    G_norm2::Vector{Float64}
    Δorb::Matrix{SVector{2,Float64}}
    G_orb_phase::Array{ComplexF64,3}
end
@inline (v::LongRangePart2D)(::Val{+}, A, k::ReducedCoordinates, args...; kwargs...) =
    v(Val(+), A, k[v.xyindex], args...; kwargs...)
function (v::LongRangePart2D)(::Val{+}, A, k::SVector{2,<:Real}; isΓ::Bool=false, nearΓ::Bool=isΓ)
    if nearΓ
        # the remain terms at k=0, only omit 1/|k|.
        for j in 2:v.norb, i in 1:(j-1)
            φᵢⱼ = v.φ[i, j]
            k_remain = φᵢⱼ(Val(:taylor_0))
            A[i, j] += k_remain
            A[j, i] += conj(k_remain)
        end
        for i in 1:v.norb
            φᵢᵢ = v.φ[i, i]
            k_remain = φᵢᵢ(Val(:taylor_0))
            A[i, i] += k_remain
        end
        return _LongRangePart2D_wok!(A, v, normalize_kdirect(k), Val(isΓ))
    else
        return _LongRangePart2D_wik!(A, v, normalize_kdirect(k))
    end
end
function (v::LongRangePart2D)(::Val{+}, A, k::SVector{2,<:Real}, nkorkgrid; isΓ::Bool=false, nearΓ::Bool=isΓ)
    if nearΓ
        _LongRangePart2D_head!(A, v, nkorkgrid)
        return _LongRangePart2D_wok!(A, v, normalize_kdirect(k), Val(isΓ))
    else
        return _LongRangePart2D_wik!(A, v, normalize_kdirect(k))
    end
end
(v::LongRangePart2D)(::Val{:G0scale}, A, k::ReducedCoordinates) = v(Val(:G0scale), A, k[v.xyindex])
function (v::LongRangePart2D)(::Val{:G0scale}, A, k::SVector{2,<:Real})
    for j in 2:v.norb, i in 1:(j-1)
        A[i, j] = v.φ[i, j](Val(:kscale), k) * cispi(2 * (k ⋅ v.Δorb[i, j]))
        A[j, i] = conj(A[i, j])
    end
    for i in 1:v.norb
        A[i, i] = v.φ[i, i](Val(:kscale), k)
    end
    return A
end
function _LongRangePart2D_wik!(A, v, k)

    k_car = v.rlattice2D * k
    kG_norm2 = map(v.G_car) do G
        kG = k_car + G
        kG ⋅ kG
    end
    k_norm2 = k_car ⋅ k_car

    φkG = Vector{Float64}(undef, v.nGminus1)
    for j in 2:v.norb, i in 1:(j-1)
        φᵢⱼ = v.φ[i, j]
        map(enumerate(kG_norm2)) do (iG, kG2)
            φkG[iG] = φᵢⱼ(Val(:k²), kG2; atol=0.0)
        end
        vlr = cispi(2 * (k ⋅ v.Δorb[i, j])) * (sum(iG -> φkG[iG] * v.G_orb_phase[iG, i, j], Base.OneTo(v.nGminus1)) +
            φᵢⱼ(Val(:k²), k_norm2; atol=0.0))
        A[i, j] += vlr
        A[j, i] += conj(vlr)
    end
    for i in 1:v.norb
        φᵢᵢ = v.φ[i, i]
        map(enumerate(kG_norm2)) do (iG, kG2)
            φkG[iG] = φᵢᵢ(Val(:k²), kG2; atol=0.0)
        end
        vlr = sum(φkG) + φᵢᵢ(Val(:k²), k_norm2; atol=0.0)
        A[i, i] += vlr
    end

    return A
end
function _LongRangePart2D_wok!(A, v, k, ::Val{false})

    # omit the term of k, but k isn't zero.
    k_car = v.rlattice2D * k
    kG_norm2 = map(v.G_car) do G
        kG = k_car + G
        kG ⋅ kG
    end

    φkG = Vector{Float64}(undef, v.nGminus1)
    for j in 2:v.norb, i in 1:(j-1)
        φᵢⱼ = v.φ[i, j]
        map(enumerate(kG_norm2)) do (iG, kG2)
            φkG[iG] = φᵢⱼ(Val(:k²), kG2; atol=0.0)
        end
        vlr = cispi(2 * (k ⋅ v.Δorb[i, j])) * sum(iG -> φkG[iG] * v.G_orb_phase[iG, i, j], Base.OneTo(v.nGminus1))
        A[i, j] += vlr
        A[j, i] += conj(vlr)
    end
    for i in 1:v.norb
        φᵢᵢ = v.φ[i, i]
        map(enumerate(kG_norm2)) do (iG, kG2)
            φkG[iG] = φᵢᵢ(Val(:k²), kG2; atol=0.0)
        end
        vlr = sum(φkG)
        A[i, i] += vlr
    end

    return A
end
function _LongRangePart2D_wok!(A, v, k, ::Val{true})

    # omit the term of k, and k is zero.

    φG = Vector{Float64}(undef, v.nGminus1)
    for j in 2:v.norb, i in 1:(j-1)
        φᵢⱼ = v.φ[i, j]
        map(enumerate(v.G_norm2)) do (iG, G2)
            φG[iG] = φᵢⱼ(Val(:k²), G2; atol=0.0)
        end
        vlr = sum(iG -> φG[iG] * real(v.G_orb_phase[iG, i, j]), Base.OneTo(v.nGminus1))
        A[i, j] += vlr
        A[j, i] += vlr
    end
    for i in 1:v.norb
        φᵢᵢ = v.φ[i, i]
        map(enumerate(v.G_norm2)) do (iG, G2)
            φG[iG] = φᵢᵢ(Val(:k²), G2; atol=0.0)
        end
        vlr = sum(φG)
        A[i, i] += vlr
    end

    return A
end
function _LongRangePart2D_head!(A, v, nk::Integer)
    signbit(nk) && error("nk should be positive")

    for j in 2:v.norb, i in 1:(j-1)
        φᵢⱼ = v.φ[i, j]
        head = φᵢⱼ(Val(:head), nk)
        A[i, j] += head
        A[j, i] += head
    end
    for i in 1:v.norb
        φᵢᵢ = v.φ[i, i]
        head = φᵢᵢ(Val(:head), nk)
        A[i, i] += head
    end

    return A
end
function _LongRangePart2D_head!(A, v, kgrid::Union{AbstractVector{<:Integer},NTuple{3,<:Integer}})
    all(signbit, kgrid) && error("nk should be positive")
    length(kgrid) == 3 || error("kgrid should have three components")
    kgrid[v.zindex] == 1 || error("kgrid should have only one point along the non-periodic direction")
    kgrid = kgrid[v.xyindex]

    minirlattice = [v.rlattice2D[1:2, 1] ./ kgrid[1] v.rlattice2D[1:2, 2] ./ kgrid[2]]
    miniFBZ = WignerSeitz(minirlattice)
    nk = prod(kgrid)

    for j in 2:v.norb, i in 1:(j-1)
        φᵢⱼ = v.φ[i, j]
        head = φᵢⱼ(Val(:head), miniFBZ, nk)
        A[i, j] += head
        A[j, i] += head
    end
    for i in 1:v.norb
        φᵢᵢ = v.φ[i, i]
        head = φᵢᵢ(Val(:head), miniFBZ, nk)
        A[i, i] += head
    end

    return A
end
