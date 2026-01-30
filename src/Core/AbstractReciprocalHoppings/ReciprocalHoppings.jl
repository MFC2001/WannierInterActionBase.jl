"""
	ReciprocalHoppings{T, U <: BaseReciprocalHoppings{T}} <: AbstractReciprocalHoppings{T}

```julia
julia> hr isa wannier90_hr
true

julia> rh = ReciprocalHoppings(hr);

julia> k isa ReducedCoordinates
true

julia> rh(k);
```
"""
struct ReciprocalHoppings{T, U <: BaseReciprocalHoppings{T}} <: AbstractReciprocalHoppings{T}
	core::U
end

@inline function (rh::ReciprocalHoppings)(k::ReducedCoordinates)
	A = Matrix{ComplexF64}(undef, numorb(rh), numorb(rh))
	return rh(A, k)
end
@inline function (rh::ReciprocalHoppings)(A::AbstractMatrix, k::ReducedCoordinates)
	size(A) == (numorb(rh), numorb(rh)) || error("Buffer size mismatch.")
	@inbounds map(CartesianIndices((numorb(rh), numorb(rh)))) do I
		A[I] = _Arh_getindex(rh, I, k)
	end
	return A
end

@inline function (rh::ReciprocalHoppings)(k::ReducedCoordinates, orblocat)
	A = Matrix{ComplexF64}(undef, numorb(rh), numorb(rh))
	return rh(A, k, orblocat)
end
@inline function (rh::ReciprocalHoppings)(A::AbstractMatrix, k::ReducedCoordinates, orblocat)
	size(A) == (numorb(rh), numorb(rh)) || error("Buffer size mismatch.")
	length(orblocat) == numorb(rh) || error("Mismatched orbital locations.")
	@inbounds map(CartesianIndices((numorb(rh), numorb(rh)))) do I
		A[I] = _Arh_getindex(rh, I, k, orblocat)
	end
	return A
end

@inline function (rh::ReciprocalHoppings)(::Val{:partial}, lattice::Lattice, k::ReducedCoordinates, orblocat)
	A = Array{ComplexF64}(undef, numorb(rh), numorb(rh), 3)
	return rh(A, Val(:partial), lattice, k, orblocat)
end
@inline function (rh::ReciprocalHoppings)(::Val{:partial}, A::AbstractArray, lattice::Lattice, k::ReducedCoordinates, orblocat)
	size(A) == (numorb(rh), numorb(rh), 3) || error("Buffer size mismatch.")
	length(orblocat) == numorb(rh) || error("Mismatched orbital locations.")
	@inbounds map(CartesianIndices((numorb(rh), numorb(rh)))) do I
		A[I, :] = _Arh_getindex(rh, I, Val(:partial), lattice, k, orblocat)
	end
	return A
end

"""
	ReciprocalHoppings(hr::wannier90_hr{T}) -> ReciprocalHoppings{T, BaseReciprocalHoppings{T}}
"""
function ReciprocalHoppings(hr::wannier90_hr{T}) where {T}
	return ReciprocalHoppings{T, BaseReciprocalHoppings{T}}(BaseReciprocalHoppings(hr))
end
function ReciprocalHoppings(Arh::AbstractReciprocalHoppings{T}) where {T}
	return ReciprocalHoppings{T, BaseReciprocalHoppings{T}}(BaseReciprocalHoppings(Arh))
end

Base.zero(::ReciprocalHoppings{T, U}, norb::Integer) where {T, U} = ReciprocalHoppings(zero(U, norb))
Base.zero(::Type{ReciprocalHoppings{T, U}}, norb::Integer) where {T, U} = ReciprocalHoppings(zero(U, norb))

function Base.union(rh₁::ReciprocalHoppings{T₁, U₁}, rh₂::ReciprocalHoppings{T₂, U₂}) where {T₁, U₁, T₂, U₂}
	return ReciprocalHoppings(union(rh₁.core, rh₂.core))
end
function Base.union(rh₁::ReciprocalHoppings{T₁, U₁}, brh₂::BaseReciprocalHoppings{T₂}) where {T₁, U₁, T₂}
	return ReciprocalHoppings(union(rh₁.core, brh₂))
end

function spinrh(rh::ReciprocalHoppings{T, U}; mode = conj) where {T, U}
	return ReciprocalHoppings(spinrh(rh.core; mode))
end
