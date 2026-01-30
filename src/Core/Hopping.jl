export Hopping, path, value, hopphase
"""
	Hopping{T <: Number}

A `Hopping` represents a hopping term in a tight-binding model,
defined by an initial orbital index `i` in unitcell `[0, 0, 0]`,
a final orbital index `j` in unitcell `R` (in reduced coordinates),
and a hopping value `t` of type `T`.

```julia
julia> hop = Hopping([1, 2, 3, 4, 5], 0.5)
t₄₅[[1, 2, 3]] = 0.5

julia> hop = Hopping(4, 5, [1, 2, 3], 0.5)
t₄₅[[1, 2, 3]] = 0.5

julia> path(hop)
3-element ReducedCoordinates{Int64} with indices SOneTo(3):
 1
 2
 3

julia> value(hop)
0.5

julia> k = ReducedCoordinates([0.1, 0.2, 0.3]);

julia> hopphase(hop, k) == exp(1im * 2π * (k ⋅path(hop)))
true
```
"""
struct Hopping{T <: Number}
	i::Int
	j::Int
	R::ReducedCoordinates{Int}
	t::T
end
"""
	Hopping(path::AbstractVector{<:Integer}, t::T) where {T <: Number} -> Hopping{T}
Creat es a `Hopping` from a vector of the form `[Rx, Ry, Rz, i, j]` and a hopping value `t`.
"""
function Hopping(path::AbstractVector{<:Integer}, t::T) where {T <: Number}
	length(path) == 5 || throw(ArgumentError("path-vector must have exactly 5 integer elements"))
	return Hopping{T}(path[4], path[5], ReducedCoordinates{Int}(path[1], path[2], path[3]), t)
end
"""
	Hopping(i::Integer, j::Integer, R::AbstractVector{<:Integer}, t::T) where {T <: Number} -> Hopping{T}
Creates a `Hopping` from integers `i`, `j`, a vector of the form `[Rx, Ry, Rz]` and a hopping value `t`.
"""
function Hopping(i::Integer, j::Integer, R::AbstractVector{<:Integer}, t::T) where {T <: Number}
	length(R) == 3 || throw(ArgumentError("R-vector must have exactly 3 integer elements"))
	return Hopping{T}(i, j, ReducedCoordinates{Int}(R), t)
end

@inline path(hop::Hopping) = hop.R
@inline value(hop::Hopping) = hop.t
function Base.show(io::IO, hop::Hopping)
	print(io, "t", subscriptnumber(hop.i), subscriptnumber(hop.j), "[", hop.R, "] = ", hop.t)
end
Base.isequal(hop₁::Hopping{T₁}, hop₂::Hopping{T₂}) where {T₁, T₂} =
	hop₁.i == hop₂.i && hop₁.j == hop₂.j && hop₁.R == hop₂.R && isequal(hop₁.t, hop₂.t)
Base.isapprox(hop₁::Hopping{T₁}, hop₂::Hopping{T₂}; atol = 1e-4) where {T₁, T₂} =
	hop₁.i == hop₂.i && hop₁.j == hop₂.j && hop₁.R == hop₂.R && isapprox(hop₁.t, hop₂.t; atol)

Base.convert(::Type{Hopping{T₁}}, hop::Hopping{T₂}) where {T₁, T₂} = Hopping{T₁}(hop.i, hop.j, hop.R, convert(T₁, hop.t))
Base.convert(::Hopping{T₁}, hop::Hopping{T₂}) where {T₁, T₂} = Hopping{T₁}(hop.i, hop.j, hop.R, convert(T₁, hop.t))
Base.promote_rule(::Type{Hopping{T₁}}, ::Type{Hopping{T₂}}) where {T₁ <: Number, T₂ <: Number} = Hopping{promote_type(T₁, T₂)}

function Base.conj(hop::Hopping{T}) where {T}
	return Hopping{T}(hop.j, hop.i, -hop.R, conj(hop.t))
end
function Base.length(hop::Hopping{T}, lattice) where {T}
	return norm(lattice * path(hop))
end
function Base.length(hop::Hopping{T}, lattice, orblocation) where {T}
	return norm(lattice * (path(hop) + orblocation[hop.j] - orblocation[hop.i]))
end
function Base.similar(hop::Hopping{T₁}, new_value::T₂) where {T₁, T₂ <: Number}
	return Hopping{T₂}(hop.i, hop.j, hop.R, new_value)
end

@inline function hopphase(hop::Hopping{T}, k::ReducedCoordinates) where {T}
	return cis(2π * (k ⋅ hop.R))
end

