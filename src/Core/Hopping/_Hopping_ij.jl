struct _Hopping_ij{S<:Integer,T<:Number}
    R::ReducedCoordinates{S}
    t::T
end
function _Hopping_ij(R::AbstractVector{S}, t::T) where {S<:Integer,T<:Number}
    length(R) == 3 || throw(ArgumentError("R-vector must have exactly 3 integer elements"))
    return _Hopping_ij{S,T}(ReducedCoordinates{S}(R), t)
end

@inline path(hop::_Hopping_ij) = hop.R
@inline value(hop::_Hopping_ij) = hop.t
@inline hopphase(hop::_Hopping_ij, k::ReducedCoordinates) = cispi(2 * (k ⋅ hop.R))

function Base.show(io::IO, hop::_Hopping_ij)
    print(io, "t", "[", hop.R, "] = ", hop.t)
end
Base.isequal(hop₁::_Hopping_ij{S₁,T₁}, hop₂::_Hopping_ij{S₂,T₂}) where {S₁,T₁,S₂,T₂} =
    hop₁.R == hop₂.R && isequal(hop₁.t, hop₂.t)
Base.hash(hop::_Hopping_ij{S,T}, h::UInt) where {S,T} = hash((hop.R, hop.t), h)
Base.isapprox(hop₁::_Hopping_ij{S₁,T₁}, hop₂::_Hopping_ij{S₂,T₂}; kwargs...) where {S₁,T₁,S₂,T₂} =
    hop₁.R == hop₂.R && isapprox(hop₁.t, hop₂.t; kwargs...)

Base.convert(::Type{_Hopping_ij{S₁,T₁}}, hop::_Hopping_ij{S₂,T₂}) where {S₁,T₁,S₂,T₂} =
    _Hopping_ij{S₁,T₁}(hop.R, convert(T₁, hop.t))
Base.promote_rule(::Type{_Hopping_ij{S₁,T₁}}, ::Type{_Hopping_ij{S₂,T₂}}) where {S₁,T₁,S₂,T₂} =
    _Hopping_ij{promote_type(S₁, S₂),promote_type(T₁, T₂)}

Base.iterate(hop::_Hopping_ij) = (hop, nothing)
Base.iterate(::_Hopping_ij, ::Any) = nothing
Base.length(::_Hopping_ij) = 1
Base.isempty(hop::_Hopping_ij) = false

Base.real(hop::_Hopping_ij{S,T}) where {S,T<:Real} = hop
Base.real(hop::_Hopping_ij{S,T}) where {S,T<:Complex} = similar(hop, real(hop.t))

Base.isreal(hop::_Hopping_ij{S,T}) where {S,T<:Real} = true
Base.isreal(hop::_Hopping_ij{S,T}) where {S,T<:Complex} = isreal(hop.t)

Base.zero(hop::_Hopping_ij{S,T}) where {S,T} = similar(hop, zero(T))
Base.iszero(hop::_Hopping_ij{S,T}) where {S,T} = iszero(hop.t)

function Base.conj(hop::_Hopping_ij{S,T}) where {S,T}
    return _Hopping_ij{S,T}(-hop.R, conj(hop.t))
end
function Base.similar(hop::_Hopping_ij{S,T₁}, new_value::T₂) where {S,T₁,T₂<:Number}
    return _Hopping_ij{S,T₂}(hop.R, new_value)
end
function distance(hop::_Hopping_ij, lattice)
    return norm(lattice * path(hop))
end
@inline _hop_type(::AbstractArray{_Hopping_ij{S,T}}) where {S,T} = _Hopping_ij{S,T}
@inline _hop_type(::_Hopping_ij{S,T}) where {S,T} = _Hopping_ij{S,T}
@inline function _hop_type(hops::AbstractArray{<:_Hopping_ij})
    isempty(hops) && return _Hopping_ij{Int8,Int}
    return reduce(promote_type, _hop_type(hop) for hop in hops)
end
function Base.unique(hops::AbstractArray{<:_Hopping_ij})
    hT = _hop_type(hops)
    (S, T) = hT.parameters
    hop_dict = Dict{ReducedCoordinates{S},T}()
    for hop in hops
        hop_dict[hop.R] = get(hop_dict, hop.R, zero(T)) + convert(T, hop.t)
    end
    unique_hops = similar(hops, _Hopping_ij{S,T}, length(hop_dict))
    for (idx, (R, t)) in enumerate(hop_dict)
        unique_hops[idx] = _Hopping_ij{S,T}(R, t)
    end
    return unique_hops
end
function Base.unique!(hops::AbstractArray{<:_Hopping_ij})
    hT = _hop_type(hops)
    (S, T) = hT.parameters
    hop_dict = Dict{ReducedCoordinates{S},T}()
    for hop in hops
        hop_dict[hop.R] = get(hop_dict, hop.R, zero(T)) + convert(T, hop.t)
    end
    empty!(hops)
    sizehint!(hops, length(hop_dict))
    for (R, t) in hop_dict
        push!(hops, _Hopping_ij{S,T}(R, t))
    end
    return hops
end







# function _ijhop_vec_to_dict!(hop_dict::Dict{ReducedCoordinates{S},T}, hop::_Hopping_ij) where {S,T}
#     key = hop.R
#     hop_dict[key] = get(hop_dict, key, zero(T)) + convert(T, hop.t)
#     return hop_dict
# end
# function _ijhop_vec_to_dict!(hop_dict::Dict{ReducedCoordinates{S},T}, hops) where {S,T}
#     for hop in hops
#         key = hop.R
#         hop_dict[key] = get(hop_dict, key, zero(T)) + convert(T, hop.t)
#     end
#     return hop_dict
# end

# function _merge_ijhops_array(hops::AbstractArray{<:_Hopping_ij})
#     isempty(hops) && return similar(hops, 0)

#     hT = _hop_type(hops)
#     (S, T) = hT.parameters
#     hop_dict = Dict{ReducedCoordinates{S},T}()
#     _ijhop_vec_to_dict!(hop_dict, hops)

#     all_hops = similar(hops, _Hopping_ij{S,T}, length(hop_dict))
#     for (idx, (R, t)) in enumerate(hop_dict)
#         all_hops[idx] = _Hopping_ij{S,T}(R, t)
#     end

#     return all_hops
# end
# function _merge_ijhops_array!(hops::AbstractArray{<:_Hopping_ij})
#     isempty(hops) && return hops

#     hT = _hop_type(hops)
#     (S, T) = hT.parameters
#     hop_dict = Dict{ReducedCoordinates{S},T}()
#     _ijhop_vec_to_dict!(hop_dict, hops)

#     empty!(hops)
#     sizehint!(hops, length(hop_dict))
#     for (R, t) in hop_dict
#         push!(hops, _Hopping_ij{S,T}(R, t))
#     end

#     return hops
# end



# function _merge_ijhops(hops::AbstractArray{<:_Hopping_ij}, arrays::Union{AbstractArray{<:_Hopping_ij},_Hopping_ij}...)
#     arrays_flat = Iterators.flatten(arrays)
#     isempty(arrays_flat) && return _merge_ijhops_array(hops)

#     hT = reduce(promote_type, _hop_type(arr) for arr in arrays; init=_hop_type(hops))
#     (S, T) = hT.parameters
#     hop_dict = Dict{ReducedCoordinates{S},T}()
#     _ijhop_vec_to_dict!(hop_dict, hops)
#     _ijhop_vec_to_dict!(hop_dict, arrays_flat)

#     all_hops = similar(hops, _Hopping_ij{S,T}, length(hop_dict))
#     for (idx, (R, t)) in enumerate(hop_dict)
#         all_hops[idx] = _Hopping_ij{S,T}(R, t)
#     end

#     return all_hops
# end
# function merge_ijhops(hop::_Hopping_ij, arrays::Union{AbstractArray{<:_Hopping_ij},_Hopping_ij}...)
#     arrays_flat = Iterators.flatten(arrays)
#     isempty(arrays_flat) && return [hop]

#     hT = reduce(promote_type, _hop_type(arr) for arr in arrays; init=_hop_type(hop))
#     (S, T) = hT.parameters
#     hop_dict = Dict{ReducedCoordinates{S},T}()
#     _ijhop_vec_to_dict!(hop_dict, hop)
#     _ijhop_vec_to_dict!(hop_dict, arrays_flat)

#     all_hops = Vector{_Hopping_ij{S,T}}(undef, length(hop_dict))
#     for (idx, (R, t)) in enumerate(hop_dict)
#         all_hops[idx] = _Hopping_ij{S,T}(R, t)
#     end

#     return all_hops
# end
# function merge_ijhops!(hops::AbstractVector{<:_Hopping_ij}, arrays::Union{AbstractArray{<:_Hopping_ij},_Hopping_ij}...)
#     arrays_flat = Iterators.flatten(arrays)
#     isempty(arrays_flat) && return _merge_ijhops_array!(hops)

#     hT = reduce(promote_type, _hop_type(arr) for arr in arrays; init=_hop_type(hops))
#     hT <: eltype(hops) || error("Incompatible types for `merge!`")
#     (S, T) = hT.parameters
#     hop_dict = Dict{ReducedCoordinates{S},T}()
#     _ijhop_vec_to_dict!(hop_dict, hops)
#     _ijhop_vec_to_dict!(hop_dict, arrays_flat)

#     empty!(hops)
#     sizehint!(hops, length(hop_dict))
#     for (R, t) in hop_dict
#         push!(hops, _Hopping_ij{S,T}(R, t))
#     end

#     return hops
# end
