export AbstractReciprocalHoppings, ReciprocalHoppings, HermitianReciprocalHoppings

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

@inline numorb(Arh::AbstractReciprocalHoppings) = numorb(Arh.core)
@inline Nhop(Arh::AbstractReciprocalHoppings) = Nhop(Arh.core)
@inline hops(Arh::AbstractReciprocalHoppings) = hops(Arh.core)
@inline BaseReciprocalHoppings(Arh::AbstractReciprocalHoppings) = BaseReciprocalHoppings(Arh.core)
wannier90_hr(Arh::AbstractReciprocalHoppings) = wannier90_hr(Arh.core)

@inline function (Arh::AbstractReciprocalHoppings)(sym::Symbol, paras...)
	return Arh(Val(sym), paras...)
end
@inline _Arh_getindex(Arh::AbstractReciprocalHoppings, paras...) = _Arh_getindex(Arh.core, paras...)

function Base.union(rh₁::AbstractReciprocalHoppings{T₁}, rh₂::AbstractReciprocalHoppings{T₂}) where {T₁, T₂}
	return ReciprocalHoppings(union(BaseReciprocalHoppings(rh₁), BaseReciprocalHoppings(rh₂)))
end
