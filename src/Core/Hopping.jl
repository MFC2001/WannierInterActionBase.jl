export Hopping, path, value, hopphase, distance
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
@inline hopphase(hop::Hopping, k::ReducedCoordinates) = cispi(2 * (k ⋅ hop.R))

function Base.show(io::IO, hop::Hopping)
	print(io, "t", subscriptnumber(hop.i), subscriptnumber(hop.j), "[", hop.R, "] = ", hop.t)
end
Base.isequal(hop₁::Hopping{T₁}, hop₂::Hopping{T₂}) where {T₁, T₂} =
	hop₁.i == hop₂.i && hop₁.j == hop₂.j && hop₁.R == hop₂.R && isequal(hop₁.t, hop₂.t)
Base.hash(hop::Hopping{T}, h::UInt) where {T} = hash((hop.i, hop.i, hop.R, hop.t), h)
Base.isapprox(hop₁::Hopping{T₁}, hop₂::Hopping{T₂}; atol = 1e-4) where {T₁, T₂} =
	hop₁.i == hop₂.i && hop₁.j == hop₂.j && hop₁.R == hop₂.R && isapprox(hop₁.t, hop₂.t; atol)

Base.convert(::Type{Hopping{T₁}}, hop::Hopping{T₂}) where {T₁, T₂} = Hopping{T₁}(hop.i, hop.j, hop.R, convert(T₁, hop.t))
Base.promote_rule(::Type{Hopping{T₁}}, ::Type{Hopping{T₂}}) where {T₁ <: Number, T₂ <: Number} = Hopping{promote_type(T₁, T₂)}

Base.iterate(hop::Hopping) = (hop, nothing)
Base.iterate(::Hopping, ::Any) = nothing
Base.length(::Hopping) = 1
Base.isempty(hop::Hopping) = false

Base.real(hop::Hopping{T}) where {T <: Real} = hop
Base.real(hop::Hopping{T}) where {T <: Complex} = similar(hop, real(hop.t))

Base.isreal(hop::Hopping{T}) where {T <: Real} = true
Base.isreal(hop::Hopping{T}) where {T <: Complex} = isreal(hop.t)

Base.zero(hop::Hopping{T}) where {T} = similar(hop, zero(T))
Base.iszero(hop::Hopping{T}) where {T} = iszero(hop.t)

function Base.conj(hop::Hopping{T}) where {T}
	return Hopping{T}(hop.j, hop.i, -hop.R, conj(hop.t))
end
function Base.similar(hop::Hopping{T₁}, new_value::T₂) where {T₁, T₂ <: Number}
	return Hopping{T₂}(hop.i, hop.j, hop.R, new_value)
end
function distance(hop::Hopping{T}, lattice) where {T}
	return norm(lattice * path(hop))
end
function distance(hop::Hopping{T}, lattice, orblocation) where {T}
	return norm(lattice * (path(hop) + orblocation[hop.j] - orblocation[hop.i]))
end

@inline _hop_amplitude_type(::AbstractArray{Hopping{T}}) where {T} = T
@inline _hop_amplitude_type(::Hopping{T}) where {T} = T
@inline function _hop_amplitude_type(hops::AbstractArray{Hopping})
	isempty(hops) && return Int
	return reduce(promote_type, _hop_amplitude_type(hop) for hop in hops)
end
function _hop_vec_to_dict!(hop_dict::Dict{Tuple{Int, Int, ReducedCoordinates{Int}}, T}, hop::Hopping) where {T}
	key = (hop.i, hop.j, hop.R)
	hop_dict[key] = get(hop_dict, key, zero(T)) + convert(T, hop.t)
	return hop_dict
end
function _hop_vec_to_dict!(hop_dict::Dict{Tuple{Int, Int, ReducedCoordinates{Int}}, T}, hops) where {T}
	for hop in hops
		key = (hop.i, hop.j, hop.R)
		hop_dict[key] = get(hop_dict, key, zero(T)) + convert(T, hop.t)
	end
	return hop_dict
end
function _merge_hops_array(hops::AbstractArray{<:Hopping})
	isempty(hops) && return similar(hops, 0)

	T = _hop_amplitude_type(hops)
	hop_dict = Dict{Tuple{Int, Int, ReducedCoordinates{Int}}, T}()
	_hop_vec_to_dict!(hop_dict, hops)

	all_hops = similar(hops, Hopping{T}, length(hop_dict))
	for (idx, (key, t)) in enumerate(hop_dict)
		i, j, R = key
		all_hops[idx] = Hopping{T}(i, j, R, t)
	end

	return all_hops
end
function _merge_hops_array!(hops::AbstractVector{<:Hopping})
	isempty(hops) && return hops

	T = _hop_amplitude_type(hops)
	hop_dict = Dict{Tuple{Int, Int, ReducedCoordinates{Int}}, T}()
	_hop_vec_to_dict!(hop_dict, hops)

	empty!(hops)
	sizehint!(hops, length(hop_dict))
	for (key, t) in hop_dict
		i, j, R = key
		push!(hops, Hopping{T}(i, j, R, t))
	end

	return hops
end
function Base.merge(hops::AbstractArray{<:Hopping}, arrays::Union{AbstractArray{<:Hopping}, Hopping}...)
	arrays_flat = Iterators.flatten(arrays)
	isempty(arrays_flat) && return _merge_hops_array(hops)

	T = reduce(promote_type, _hop_amplitude_type(arr) for arr in arrays; init = _hop_amplitude_type(hops))

	hop_dict = Dict{Tuple{Int, Int, ReducedCoordinates{Int}}, T}()
	_hop_vec_to_dict!(hop_dict, hops)
	_hop_vec_to_dict!(hop_dict, arrays_flat)

	all_hops = similar(hops, Hopping{T}, length(hop_dict))
	for (idx, (key, t)) in enumerate(hop_dict)
		i, j, R = key
		all_hops[idx] = Hopping{T}(i, j, R, t)
	end

	return all_hops
end
function Base.merge(hop::Hopping, arrays::Union{AbstractArray{<:Hopping}, Hopping}...)
	arrays_flat = Iterators.flatten(arrays)
	isempty(arrays_flat) && return [hop]

	T = reduce(promote_type, _hop_amplitude_type(arr) for arr in arrays; init = _hop_amplitude_type(hop))

	hop_dict = Dict{Tuple{Int, Int, ReducedCoordinates{Int}}, T}()
	_hop_vec_to_dict!(hop_dict, hop)
	_hop_vec_to_dict!(hop_dict, arrays_flat)

	all_hops = Vector{Hopping{T}}(undef, length(hop_dict))
	for (idx, (key, t)) in enumerate(hop_dict)
		i, j, R = key
		all_hops[idx] = Hopping{T}(i, j, R, t)
	end

	return all_hops
end
function Base.merge!(hops::AbstractVector{<:Hopping}, arrays::Union{AbstractArray{<:Hopping}, Hopping}...)
	arrays_flat = Iterators.flatten(arrays)
	isempty(arrays_flat) && return _merge_hops_array!(hops)

	T = reduce(promote_type, _hop_amplitude_type(arr) for arr in arrays; init = _hop_amplitude_type(hops))
	Hopping{T} <: eltype(hops) || error("Incompatible types for `merge!`")

	hop_dict = Dict{Tuple{Int, Int, ReducedCoordinates{Int}}, T}()
	_hop_vec_to_dict!(hop_dict, hops)
	_hop_vec_to_dict!(hop_dict, arrays_flat)

	empty!(hops)
	sizehint!(hops, length(hop_dict))
	for (key, t) in hop_dict
		i, j, R = key
		push!(hops, Hopping{T}(i, j, R, t))
	end

	return hops
end
function _ijhop_vec_to_dict!(hop_dict::Dict{ReducedCoordinates{Int}, T}, hop::Hopping) where {T}
	key = hop.R
	hop_dict[key] = get(hop_dict, key, zero(T)) + convert(T, hop.t)
	return hop_dict
end
function _ijhop_vec_to_dict!(hop_dict::Dict{ReducedCoordinates{Int}, T}, hops) where {T}
	for hop in hops
		key = hop.R
		hop_dict[key] = get(hop_dict, key, zero(T)) + convert(T, hop.t)
	end
	return hop_dict
end
function _merge_ijhops_array(hops::AbstractArray{<:Hopping})
	isempty(hops) && return similar(hops, 0)

	T = _hop_amplitude_type(hops)
	hop_dict = Dict{ReducedCoordinates{Int}, T}()
	_ijhop_vec_to_dict!(hop_dict, hops)

	i, j = first(hops).i, first(hops).j
	all_hops = similar(hops, Hopping{T}, length(hop_dict))
	for (idx, (R, t)) in enumerate(hop_dict)
		all_hops[idx] = Hopping{T}(i, j, R, t)
	end

	return all_hops
end
function _merge_ijhops_array!(hops::AbstractVector{<:Hopping})
	isempty(hops) && return hops

	T = _hop_amplitude_type(hops)
	hop_dict = Dict{ReducedCoordinates{Int}, T}()
	_ijhop_vec_to_dict!(hop_dict, hops)

	i, j = first(hops).i, first(hops).j
	empty!(hops)
	sizehint!(hops, length(hop_dict))
	for (R, t) in hop_dict
		push!(hops, Hopping{T}(i, j, R, t))
	end

	return hops
end
function merge_ijhops(hops::AbstractArray{<:Hopping}, arrays::Union{AbstractArray{<:Hopping}, Hopping}...)
	arrays_flat = Iterators.flatten(arrays)
	isempty(arrays_flat) && return _merge_ijhops_array(hops)

	T = reduce(promote_type, _hop_amplitude_type(arr) for arr in arrays; init = _hop_amplitude_type(hops))

	hop_dict = Dict{ReducedCoordinates{Int}, T}()
	_ijhop_vec_to_dict!(hop_dict, hops)
	_ijhop_vec_to_dict!(hop_dict, arrays_flat)

	i, j = first(arrays_flat).i, first(arrays_flat).j
	all_hops = similar(hops, Hopping{T}, length(hop_dict))
	for (idx, (R, t)) in enumerate(hop_dict)
		all_hops[idx] = Hopping{T}(i, j, R, t)
	end

	return all_hops
end
function merge_ijhops(hop::Hopping, arrays::Union{AbstractArray{<:Hopping}, Hopping}...)
	arrays_flat = Iterators.flatten(arrays)
	isempty(arrays_flat) && return [hop]

	T = reduce(promote_type, _hop_amplitude_type(arr) for arr in arrays; init = _hop_amplitude_type(hop))

	hop_dict = Dict{ReducedCoordinates{Int}, T}()
	_ijhop_vec_to_dict!(hop_dict, hop)
	_ijhop_vec_to_dict!(hop_dict, arrays_flat)

	i, j = first(arrays_flat).i, first(arrays_flat).j
	all_hops = Vector{Hopping{T}}(undef, length(hop_dict))
	for (idx, (R, t)) in enumerate(hop_dict)
		all_hops[idx] = Hopping{T}(i, j, R, t)
	end

	return all_hops
end
function merge_ijhops!(hops::AbstractVector{<:Hopping}, arrays::Union{AbstractArray{<:Hopping}, Hopping}...)
	arrays_flat = Iterators.flatten(arrays)
	isempty(arrays_flat) && return _merge_ijhops_array!(hops)

	T = reduce(promote_type, _hop_amplitude_type(arr) for arr in arrays; init = _hop_amplitude_type(hops))
	Hopping{T} <: eltype(hops) || error("Incompatible types for `merge!`")

	hop_dict = Dict{ReducedCoordinates{Int}, T}()
	_ijhop_vec_to_dict!(hop_dict, hops)
	_ijhop_vec_to_dict!(hop_dict, arrays_flat)

	i, j = first(arrays_flat).i, first(arrays_flat).j
	empty!(hops)
	sizehint!(hops, length(hop_dict))
	for (R, t) in hop_dict
		push!(hops, Hopping{T}(i, j, R, t))
	end

	return hops
end
