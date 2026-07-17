"""
	ReciprocalHoppings{T, U <: BaseReciprocalHoppings{T}} <: AbstractReciprocalHoppings{T}

```julia
julia> hr isa HR
true

julia> rh = ReciprocalHoppings(hr);

julia> k isa ReducedCoordinates
true

julia> rh(k);
```
"""
struct ReciprocalHoppings{T, U <: BaseReciprocalHoppings{T}} <: AbstractReciprocalHoppings{T}
	core::U
	ij::Vector{Tuple{Int, Int}}
end

@inline function (rh::ReciprocalHoppings)(k::ReducedCoordinates)
	A = Matrix{ComplexF64}(undef, numorb(rh), numorb(rh))
	_Arh_calc!(A, rh.core, rh.ij, k)
	return A
end
@inline function (rh::ReciprocalHoppings)(A::AbstractMatrix, k::ReducedCoordinates)
	size(A) == (numorb(rh), numorb(rh)) || error("Buffer size mismatch.")
	_Arh_calc!(A, rh.core, rh.ij, k)
	return A
end

@inline function (rh::ReciprocalHoppings)(k::ReducedCoordinates, orblocat::AbstractVector{<:ReducedCoordinates})
	A = Matrix{ComplexF64}(undef, numorb(rh), numorb(rh))
	length(orblocat) == numorb(rh) || error("Mismatched orbital locations.")
	_Arh_calc!(A, rh.core, rh.ij, k, orblocat)
	return A
end
@inline function (rh::ReciprocalHoppings)(A::AbstractMatrix, k::ReducedCoordinates, orblocat::AbstractVector{<:ReducedCoordinates})
	size(A) == (numorb(rh), numorb(rh)) || error("Buffer size mismatch.")
	length(orblocat) == numorb(rh) || error("Mismatched orbital locations.")
	_Arh_calc!(A, rh.core, rh.ij, k, orblocat)
	return A
end

@inline function (rh::ReciprocalHoppings)(::Val{:partial}, lattice::Lattice, k::ReducedCoordinates, orblocat::AbstractVector{<:ReducedCoordinates})
	A = Array{ComplexF64}(undef, numorb(rh), numorb(rh), 3)
	length(orblocat) == numorb(rh) || error("Mismatched orbital locations.")
	_Arh_calc!(A, rh.core, rh.ij, Val(:partial), k, lattice, orblocat)
	return A
end
@inline function (rh::ReciprocalHoppings)(::Val{:partial}, A::AbstractArray, lattice::Lattice, k::ReducedCoordinates, orblocat::AbstractVector{<:ReducedCoordinates})
	size(A) == (numorb(rh), numorb(rh), 3) || error("Buffer size mismatch.")
	length(orblocat) == numorb(rh) || error("Mismatched orbital locations.")
	_Arh_calc!(A, rh.core, rh.ij, Val(:partial), k, lattice, orblocat)
	return A
end

function Base.show(io::IO, rh::ReciprocalHoppings)
	print(io, "ReciprocalHoppings with $(sum(Nhop(rh))) hoppings and $(numorb(rh)) orbitals.")
end
function Base.filter(F, rh::ReciprocalHoppings)
	brh = filter(F, rh.core)
	return ReciprocalHoppings(brh)
end
function prune(rh::ReciprocalHoppings, cutoff::Real)
	brh = prune(rh.core, cutoff)
	return ReciprocalHoppings(brh)
end
function Base.real(rh::ReciprocalHoppings{T, U}) where {T <: Complex, U}
	brh = real(rh.core)
	return ReciprocalHoppings(brh)
end
Base.zero(::Type{ReciprocalHoppings{T, U}}, norb::Integer) where {T, U} =
	ReciprocalHoppings(zero(U, norb))
Base.convert(::Type{ReciprocalHoppings{T₁, U₁}}, rh::ReciprocalHoppings{T₂, U₂}) where {T₁, U₁, T₂, U₂} =
	ReciprocalHoppings(convert(U₁, rh.core))
"""
	translate(rh::ReciprocalHoppings, centres::Pair{<:Integer, <:AbstractVector{<:Integer}}...)

You only need to input the wannier centres R which isn't in the [0, 0, 0] cell.

See also `translate!`.
"""
function translate(rh::ReciprocalHoppings, aim_centres::Pair{<:Integer, <:AbstractVector{<:Integer}}...)
	brh = translate(rh.core, aim_centres)
	return ReciprocalHoppings(brh)
end
Base.merge(rh::ReciprocalHoppings, Arh::AbstractReciprocalHoppings) =
	ReciprocalHoppings(merge(rh.core, BaseReciprocalHoppings(Arh)))
Base.merge(rh::ReciprocalHoppings, hr::HR) =
	ReciprocalHoppings(merge(rh.core, hr))
Base.merge(rh::ReciprocalHoppings, arrays::Union{AbstractArray{<:Hopping}, Hopping}...) =
	ReciprocalHoppings(merge(rh.core, arrays...))

function spinrh(rh::ReciprocalHoppings{T, U}; mode = conj) where {T, U}
	return ReciprocalHoppings(spinrh(rh.core; mode))
end

"""
	ReciprocalHoppings(hr::HR{T}) -> ReciprocalHoppings{T, BaseReciprocalHoppings{T}}
"""
function ReciprocalHoppings(hr::HR{T}) where {T}
	brh = BaseReciprocalHoppings(hr)
	ij = _matrix_index(numorb(brh))
	return ReciprocalHoppings{T, BaseReciprocalHoppings{T}}(brh, ij)
end
function ReciprocalHoppings(int::WanIntijklR{T}) where {T}
	brh = BaseReciprocalHoppings(int)
	ij = _matrix_index(numorb(brh))
	return ReciprocalHoppings{T, BaseReciprocalHoppings{T}}(brh, ij)
end
function ReciprocalHoppings(Arh::AbstractReciprocalHoppings{T}) where {T}
	brh = BaseReciprocalHoppings(Arh)
	ij = _matrix_index(numorb(brh))
	return ReciprocalHoppings{T, BaseReciprocalHoppings{T}}(brh, ij)
end
_matrix_index(norb) = reshape([(i, j) for i in 1:norb, j in 1:norb], :)
