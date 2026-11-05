export AbstractReciprocalHoppings, ReciprocalHoppings, HermitianReciprocalHoppings
export numorb, Nhop, hops, prune, prune!, shift_fermi_energy!, translate, translate!, spinorize

"""
	AbstractReciprocalHoppings{T <: Number}

Reorganized storage of the hopping terms, grouped by orbital pair, that facilitates the
evaluation of the reciprocal-space Hamiltonian. Use the concrete
[`ReciprocalHoppings`](@ref) or [`HermitianReciprocalHoppings`](@ref) types.
"""
abstract type AbstractReciprocalHoppings{T<:Number} end

#TODO
# struct BaseReciprocalHoppings{T} <: AbstractReciprocalHoppings{T} end
# struct ReciprocalHoppings{T <: Number, U :: BaseReciprocalHoppings{T}} <: AbstractReciprocalHoppings{T} end
# struct SymmetricReciprocalHoppings{T <: Number, U <: AbstractReciprocalHoppings{T}} <: AbstractReciprocalHoppings{T} end
# struct ΓReciprocalHoppings{T <: Number, U <: AbstractReciprocalHoppings{T}} <: AbstractReciprocalHoppings{T} end
# struct HermitianReciprocalHoppings{T <: Number, U <: AbstractReciprocalHoppings{T}} <: AbstractReciprocalHoppings{T} end
# function norb(arh::AbstractReciprocalHoppings) norb(arh.core) end
# function Nhop(arh::AbstractReciprocalHoppings) Nhop(arh.core) end
# function hops(arh::AbstractReciprocalHoppings) hops(arh.core) end


include("./BaseReciprocalHoppings.jl")
include("./ReciprocalHoppings.jl")
include("./HermitianReciprocalHoppings.jl")

# include("./GammaReciprocalHoppings.jl")
# include("./SymmetricReciprocalHoppings.jl")

@inline BaseReciprocalHoppings(Arh::AbstractReciprocalHoppings) = BaseReciprocalHoppings(Arh.core)
@inline numorb(Arh::AbstractReciprocalHoppings) = numorb(BaseReciprocalHoppings(Arh))
"""
	Nhop(arh::AbstractReciprocalHoppings) -> Matrix{Int}

Return the matrix `Nhop[i, j]` with the number of hopping terms stored between orbitals
`i` and `j`.
"""
@inline Nhop(Arh::AbstractReciprocalHoppings) = Nhop(BaseReciprocalHoppings(Arh))
"""
	hops(arh::AbstractReciprocalHoppings) -> Matrix{Vector}

Return the matrix whose `[i, j]` entry is the vector of hopping terms between orbitals
`i` and `j`.
"""
@inline hops(Arh::AbstractReciprocalHoppings) = hops(BaseReciprocalHoppings(Arh))
Base.eltype(::Type{AbstractReciprocalHoppings{T}}) where {T} = T
@inline Base.real(Arh::AbstractReciprocalHoppings{T}) where {T<:Real} = Arh
@inline Base.isreal(Arh::AbstractReciprocalHoppings{T}) where {T<:Real} = true
@inline Base.isreal(Arh::AbstractReciprocalHoppings{T}) where {T<:Complex} = isreal(Arh.core)
@inline Base.iszero(Arh::AbstractReciprocalHoppings) = iszero(Arh.core)
@inline function Base.filter!(F, Arh::AbstractReciprocalHoppings)
    filter!(F, Arh.core)
    return Arh
end
"""
	prune!(arh::AbstractReciprocalHoppings, cutoff::Real) -> AbstractReciprocalHoppings

Remove, in place, the hopping terms with `abs(t) < cutoff` (for complex values
`abs2(t) < cutoff^2`) and return `arh`.
"""
prune!(Arh::AbstractReciprocalHoppings, cutoff::Real) = prune!(Arh.core, cutoff)
function Base.merge!(Arh₁::AbstractReciprocalHoppings, Arh₂::AbstractReciprocalHoppings)
    merge!(Arh₁.core, Arh₂.core)
    return Arh₁
end
function Base.merge!(Arh::AbstractReciprocalHoppings, hr::AbstractHR)
    merge!(Arh.core, hr)
    return Arh
end
function Base.merge!(Arh::AbstractReciprocalHoppings, arrays::Union{AbstractArray{<:Hopping},Hopping}...)
    merge!(Arh.core, arrays...)
    return Arh
end
# HRFull(Arh::AbstractReciprocalHoppings) = HRFull(Arh.core)

@inline (Arh::AbstractReciprocalHoppings)(sym::Symbol, paras...) = Arh(Val(sym), paras...)
@inline function (Arh::AbstractReciprocalHoppings)(k::AbstractVector{<:Real})
    A = Matrix{ComplexF64}(undef, numorb(Arh), numorb(Arh))
    return Arh(A, ReducedCoordinates(k))
end
@inline function (Arh::AbstractReciprocalHoppings)(::Val{:partial}, k::AbstractVector{<:Real}, lattice::Union{Lattice,AbstractMatrix{<:Real}})
    A = Array{ComplexF64}(undef, numorb(Arh), numorb(Arh), 3)
    return Arh(A, Val(:partial), ReducedCoordinates(k), Lattice(parent(lattice)))
end
@inline function (Arh::AbstractReciprocalHoppings)(k::AbstractVector{<:Real}, orblocat::AbstractVector{<:AbstractVector{<:Real}})
    A = Matrix{ComplexF64}(undef, numorb(Arh), numorb(Arh))
    return Arh(A, ReducedCoordinates(k), ReducedCoordinates.(orblocat))
end
@inline function (Arh::AbstractReciprocalHoppings)(::Val{:partial}, k::AbstractVector{<:Real}, lattice::Union{Lattice,AbstractMatrix{<:Real}}, orblocat::AbstractVector{<:AbstractVector{<:Real}})
    A = Array{ComplexF64}(undef, numorb(Arh), numorb(Arh), 3)
    return Arh(A, Val(:partial), ReducedCoordinates(k), Lattice(parent(lattice)), ReducedCoordinates.(orblocat))
end
# """
# 	shift_fermi_energy!(hr::HR{T}, μ::Real)::HR{T} where{T}
# """
# function shift_fermi_energy!(hr::HR{T}, μ::Real)::HR{T} where {T}
#     I = findall(x -> iszero(x[1:3]) && isequal(x[4], x[5]), eachcol(hr.path))
# 	hr.value[I] .-= μ
# 	if length(I) < numorb(hr)
# 		exist_onsite = hr.path[4, I]
# 		append_onsite = Int[]
# 		for i in 1:numorb(hr)
# 			i ∈ exist_onsite && continue
# 			push!(append_onsite, i)
# 		end
# 		n_append_onsite = length(append_onsite)
# 		append!(hr)
# 	end
#     return hr
# end
"""
	translate!(rh::AbstractReciprocalHoppings{T}, centres::Pair{<:Integer, <:AbstractVector{<:Integer}}...)

You only need to input the wannier centres R which isn't in the [0, 0, 0] cell.

See also `translate`.
"""
function translate!(Arh::AbstractReciprocalHoppings{T}, aim_centres::Pair{<:Integer,<:AbstractVector{<:Integer}}...) where {T}
    norb = numorb(Arh)
    centres = [ReducedCoordinates{Int}(0, 0, 0) for _ in Base.OneTo(norb)]
    for (iw, R) in aim_centres
        centres[iw] = ReducedCoordinates{Int}(R)
    end
    for j in Base.OneTo(norb), i in Base.OneTo(norb)
        Nhop_ij = Nhop(Arh)[i, j]
        if iszero(Nhop_ij)
            continue
        end
        dR = centres[j] - centres[i]
        if iszero(dR)
            continue
        else
            hops_ij = hops(Arh)[i, j]
            for ihop in Base.OneTo(Nhop_ij)
                hop = hops_ij[ihop]
                hops_ij[ihop] = _Hopping_ij(hop.R + dR, hop.t)
            end
        end
    end
    return Arh
end
