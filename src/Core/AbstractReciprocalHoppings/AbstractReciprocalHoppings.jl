export AbstractReciprocalHoppings, ReciprocalHoppings, HermitianReciprocalHoppings
export numorb, Nhop, hops, prune, prune!, translate, translate!, spinrh

"""
	AbstractReciprocalHoppings{T <: Number}

Reorganized the storage form of the hopping terms to facilitate calculation.

	BaseReciprocalHoppings{T} <: AbstractReciprocalHoppings{T}

A concrete type of `AbstractReciprocalHoppings`. But don't use `BaseReciprocalHoppings` directly.
"""
abstract type AbstractReciprocalHoppings{T <: Number} end

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
@inline Nhop(Arh::AbstractReciprocalHoppings) = Nhop(BaseReciprocalHoppings(Arh))
@inline hops(Arh::AbstractReciprocalHoppings) = hops(BaseReciprocalHoppings(Arh))
@inline Base.real(Arh::AbstractReciprocalHoppings{T}) where {T <: Real} = Arh
@inline Base.isreal(Arh::AbstractReciprocalHoppings{T}) where {T <: Real} = true
@inline Base.isreal(Arh::AbstractReciprocalHoppings{T}) where {T <: Complex} = isreal(Arh.core)
@inline Base.iszero(Arh::AbstractReciprocalHoppings) = iszero(Arh.core)
@inline function Base.filter!(F, Arh::AbstractReciprocalHoppings)
	filter!(F, Arh.core)
	return Arh
end
prune!(Arh::AbstractReciprocalHoppings, cutoff::Real) = prune!(Arh.core, cutoff)
function Base.merge!(Arh₁::AbstractReciprocalHoppings, Arh₂::AbstractReciprocalHoppings)
	merge!(Arh₁.core, Arh₂.core)
	return Arh₁
end
function Base.merge!(Arh::AbstractReciprocalHoppings, hr::HR)
	merge!(Arh.core, hr)
	return Arh
end
function Base.merge!(Arh::AbstractReciprocalHoppings, arrays::Union{AbstractArray{<:Hopping}, Hopping}...)
	merge!(Arh.core, arrays...)
	return Arh
end
HR(Arh::AbstractReciprocalHoppings) = HR(Arh.core)

@inline (Arh::AbstractReciprocalHoppings)(sym::Symbol, paras...) = Arh(Val(sym), paras...)
@inline function (Arh::AbstractReciprocalHoppings)(k::AbstractVector{<:Real})
	A = Matrix{ComplexF64}(undef, numorb(Arh), numorb(Arh))
	return Arh(A, ReducedCoordinates(k))
end
@inline function (Arh::AbstractReciprocalHoppings)(k::AbstractVector{<:Real}, orblocat::AbstractVector{<:ReducedCoordinates})
	A = Matrix{ComplexF64}(undef, numorb(Arh), numorb(Arh))
	return Arh(A, ReducedCoordinates(k), orblocat)
end
@inline function (Arh::AbstractReciprocalHoppings)(k::AbstractVector{<:Real}, orblocat::AbstractVector{<:AbstractVector{<:Real}})
	A = Matrix{ComplexF64}(undef, numorb(Arh), numorb(Arh))
	return Arh(A, ReducedCoordinates(k), ReducedCoordinates.(orblocat))
end
@inline function (Arh::AbstractReciprocalHoppings)(::Val{:partial}, lattice::Union{Lattice, AbstractMatrix{<:Real}}, k::AbstractVector{<:Real}, orblocat::AbstractVector{<:ReducedCoordinates})
	A = Array{ComplexF64}(undef, numorb(Arh), numorb(Arh), 3)
	return Arh(A, Val(:partial), Lattice(parent(lattice)), ReducedCoordinates(k), orblocat)
end
@inline function (Arh::AbstractReciprocalHoppings)(::Val{:partial}, lattice::Union{Lattice, AbstractMatrix{<:Real}}, k::AbstractVector{<:Real}, orblocat::AbstractVector{<:AbstractVector{<:Real}})
	A = Array{ComplexF64}(undef, numorb(Arh), numorb(Arh), 3)
	return Arh(A, Val(:partial), Lattice(parent(lattice)), ReducedCoordinates(k), ReducedCoordinates.(orblocat))
end

"""
	translate!(rh::AbstractReciprocalHoppings{T}, centres::Pair{<:Integer, <:AbstractVector{<:Integer}}...)

You only need to input the wannier centres R which isn't in the [0, 0, 0] cell.

See also `translate`.
"""
function translate!(Arh::AbstractReciprocalHoppings{T}, aim_centres::Pair{<:Integer, <:AbstractVector{<:Integer}}...) where {T}
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
				hops_ij[ihop] = Hopping(i, j, hop.R + dR, hop.t)
			end
		end
	end
	return Arh
end
