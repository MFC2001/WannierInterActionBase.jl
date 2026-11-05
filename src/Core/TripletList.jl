export TripleList, vsize, linear_index, triple_index
"""
    TripleList(v1, v2, v3)

Construct a `TripleList` from three vectors `v1`, `v2` and `v3`.

The resulting object stores all combinations `(v1[i], v2[j], v3[k])` in
row‑major order (first index varies fastest). It provides fast forward
lookup by linear or triple indices, as well as reverse mapping from
values to their indices via dictionaries.

# Supported operations
- Iteration (`for x in tl`)
- `length(tl)`, `size(tl)`

# Examples
```julia
julia> tl = TripleList([1, 2], 3:5, [3, 9])
TripleList(2, 3, 2)

julia> length(tl) == 2 * 3 * 2
true

julia> tl[2]
(2, 3, 3)

julia> tl[(2, 3, 3)]   # returns the linear index of value triple (2,3,3)
2
```
"""
struct TripleList{T₁,T₂,T₃}
    values::Tuple{Vector{T₁},Vector{T₂},Vector{T₃}}
    values2idx::Tuple{Dict{T₁,Int},Dict{T₂,Int},Dict{T₃,Int}}
    data_idx::Vector{NTuple{3,Int}} # [i1, i2, i3]
    data::Vector{Tuple{T₁,T₂,T₃}} # [n1, n2, n3]
end
Base.length(tl::TripleList) = length(tl.data)
Base.size(tl::TripleList) = size(tl.data)
Base.getindex(tl::TripleList, i::Integer) = tl.data[i]
Base.getindex(tl::TripleList, I::AbstractUnitRange{<:Integer}) = tl.data[I]
Base.getindex(tl::TripleList, I::AbstractVector{<:Integer}) = tl.data[I]
Base.getindex(tl::TripleList, I::AbstractVector{Bool}) = tl.data[I]
Base.getindex(tl::TripleList, ::Colon) = tl.data[:]
Base.view(tl::TripleList, I...) = view(tl.data, I...)
function Base.iterate(tl::TripleList, i::Integer=1)
    i > length(tl) && return nothing
    return (tl[i], i + 1)
end
Base.eltype(::Type{<:TripleList{T1,T2,T3}}) where {T1,T2,T3} = Tuple{T1,T2,T3}
Base.firstindex(::TripleList) = 1
Base.lastindex(tl::TripleList) = length(tl)
Base.eachindex(tl::TripleList) = Base.OneTo(length(tl.data))
Base.keys(tl::TripleList) = Base.OneTo(length(tl.data))
Base.axes(tl::TripleList) = (Base.OneTo(length(tl.data)),)
Base.axes(tl::TripleList, d::Integer) = d == 1 ? Base.OneTo(length(tl.data)) : Base.OneTo(1)
# equal
Base.:(==)(a::TripleList, b::TripleList) = a.values == b.values
Base.hash(tl::TripleList, h::UInt) = hash(tl.values, h)
Base.isempty(tl::TripleList) = isempty(tl.data)
Base.copy(tl::TripleList) = TripleList(
    (copy(tl.values[1]), copy(tl.values[2]), copy(tl.values[3])),
    (copy(tl.values2idx[1]), copy(tl.values2idx[2]), copy(tl.values2idx[3])),
    copy(tl.data_idx), copy(tl.data))
Base.values(tl::TripleList) = tl.values
function Base.show(io::IO, tl::TripleList)
    (n1, n2, n3) = vsize(tl)
    print(io, "TripleList(", n1, ", ", n2, ", ", n3, ")")
end
"""
    vsize(tl::TripleList)
"""
vsize(tl::TripleList) = length.(tl.values)
"""
	linear_index(tl::TripleList, i::Integer, j::Integer, k::Integer) -> Int
	linear_index(tl::TripleList, value_tuple::Tuple)

Return the row-major linear index of the combination `(i, j, k)`, i.e. the index of the
value triple `value_tuple`, with the first index varying fastest.
"""
function linear_index(tl::TripleList, value_tuple::Tuple{T₁,T₂,T₃}) where {T₁,T₂,T₃}
    i1 = tl.values2idx[1][value_tuple[1]]
    i2 = tl.values2idx[2][value_tuple[2]]
    i3 = tl.values2idx[3][value_tuple[3]]
    n1 = length(tl.values[1])
    n2 = length(tl.values[2])
    return i1 + (i2-1) * n1 + (i3 - 1) * n1 * n2
end
function linear_index(tl::TripleList, i::Integer, j::Integer, k::Integer)
    n1 = length(tl.values[1])
    n2 = length(tl.values[2])
    n3 = length(tl.values[3])
    @boundscheck begin
        (1 ≤ i ≤ n1) || throw(BoundsError(tl, (i, j, k)))
        (1 ≤ j ≤ n2) || throw(BoundsError(tl, (i, j, k)))
        (1 ≤ k ≤ n3) || throw(BoundsError(tl, (i, j, k)))
    end
    return i + (j-1)*n1 + (k-1)*n1*n2
end
"""
	triple_index(tl::TripleList, idx::Integer) -> NTuple{3, Int}

Return the indices `(i, j, k)` of the `idx`-th combination stored in `tl`.
"""
@inline triple_index(tl::TripleList, idx::Integer) = tl.data_idx[idx]
# function Base.reshape(A::AbstractVector, tl::TripleList)
#     length(A) == length(tl) || throw(DimensionMismatch("length mismatch: $(length(A)) vs $(length(tl))"))
#     return reshape(A, vsize(tl))
# end
# function Base.reshape(A::AbstractArray{T,3}, tl::TripleList) where {T}
#     size(A) == vsize(tl) || throw(DimensionMismatch("size mismatch: $(size(A)) vs $(vsize(tl))"))
#     return reshape(A, :)
# end
"""
	TripleList(v1, v2, v3) -> TripleList`

Here `v1`, `v2` and `v3` should be vector.
"""
function TripleList(v1::AbstractVector{T₁}, v2::AbstractVector{T₂}, v3::AbstractVector{T₃}) where {T₁,T₂,T₃}

    v1 = collect(v1)
    v2 = collect(v2)
    v3 = collect(v3)
    values = (v1, v2, v3)

    v12idx = Dict{T₁,Int}(val => idx for (idx, val) in enumerate(v1))
    v22idx = Dict{T₂,Int}(val => idx for (idx, val) in enumerate(v2))
    v32idx = Dict{T₃,Int}(val => idx for (idx, val) in enumerate(v3))
    values2idx = (v12idx, v22idx, v32idx)

    N = length(v1) * length(v2) * length(v3)
    data_idx = Vector{NTuple{3,Int}}(undef, N)
    data = Vector{Tuple{T₁,T₂,T₃}}(undef, N)
    n = 0
    for (idx3, val3) in enumerate(v3), (idx2, val2) in enumerate(v2), (idx1, val1) in enumerate(v1)
        n += 1
        data_idx[n] = (idx1, idx2, idx3)
        data[n] = (val1, val2, val3)
    end

    return TripleList(values, values2idx, data_idx, data)
end
