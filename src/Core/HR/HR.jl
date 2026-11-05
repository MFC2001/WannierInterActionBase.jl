export AbstractHR, HRFull, HRHermitian, HR
export numorb, index, indices, prune, spinorize, translate
"""
    AbstractHR{T<:Number}

A collection of data of the form `v_{ijR}`, where `i` and `j` are integer indices
and `R` is a three‑dimensional integer vector. The numerical values `v` are of type `T`.
"""
abstract type AbstractHR{T<:Number} end
# 支持访问和修改元素，但是该操作较为低效，仅用作检查用。

include("./Full.jl")
include("./Hermitian.jl")

function Base.getproperty(hr::AbstractHR, name::Symbol)
    if name === :nR
        return length(getfield(hr, :R))
    elseif name === :nv
        return length(getfield(hr, :vals))
    else
        return getfield(hr, name)
    end
end
function _HR_Rindex(hr::AbstractHR, iR::Integer)::Int
    1 ≤ iR ≤ hr.nR || error("iR = $iR should be in [1, $(hr.nR)]")
    return iR
end
_HR_Rindex(hr::AbstractHR, R::AbstractVector{<:Integer})::Int = get(hr.R2idx, ReducedCoordinates{Int}(R), 0)
_HR_Rindex(hr::AbstractHR, R::ReducedCoordinates{Int})::Int = get(hr.R2idx, R, 0)
_HR_Rindex!(hr::AbstractHR, iR::Integer) = _HR_Rindex(hr, iR)
_HR_Rindex!(hr::AbstractHR, R::AbstractVector{<:Integer})::Int = _HR_Rindex!(hr, ReducedCoordinates{Int}(R))
# the Bool means `if_check_ij`
function _HR_getindex(hr::AbstractHR{T}, ::Val{true}, i::Integer, j::Integer, iRorR) where {T}
    1 ≤ i ≤ numorb(hr) || error("i = $i should be in [1, $(numorb(hr))]")
    1 ≤ j ≤ numorb(hr) || error("i = $j should be in [1, $(numorb(hr))]")
    iR = _HR_Rindex(hr, iRorR)
    iszero(iR) && return zero(T)
    return _HR_getindex(hr, i, j, iR)
end
function _HR_getindex(hr::AbstractHR{T}, ::Val{false}, i::Integer, j::Integer, iRorR) where {T}
    iR = _HR_Rindex(hr, iRorR)
    iszero(iR) && return zero(T)
    return _HR_getindex(hr, i, j, iR)
end
function _HR_setindex!(hr::AbstractHR{T}, ::Val{true}, v, i::Integer, j::Integer, iRorR) where {T}
    1 ≤ i ≤ numorb(hr) || error("i = $i should be in [1, $(numorb(hr))]")
    1 ≤ j ≤ numorb(hr) || error("i = $j should be in [1, $(numorb(hr))]")
    iR = _HR_Rindex!(hr, iRorR)
    return _HR_setindex!(hr, v, i, j, iR)
end
function _HR_setindex!(hr::AbstractHR{T}, ::Val{false}, v, i::Integer, j::Integer, iRorR) where {T}
    iR = _HR_Rindex!(hr, iRorR)
    return _HR_setindex!(hr, v, i, j, iR)
end
Base.getindex(hr::AbstractHR, i::Integer) = hr.vals[i]
Base.getindex(hr::AbstractHR, i::Integer, j::Integer, iRorR) = _HR_getindex(hr, Val(true), i, j, iRorR)
Base.setindex!(hr::AbstractHR{<:Real}, v::Real, i::Integer) = (hr.vals[i] = v)
Base.setindex!(hr::AbstractHR{<:Complex}, v::Number, i::Integer) = (hr.vals[i] = v)
Base.setindex!(hr::AbstractHR{<:Real}, v::Union{Tuple{<:Integer,<:Real},Real}, i, j, iRorR) = _HR_setindex!(hr, Val(true), v, i, j, iRorR)
Base.setindex!(hr::AbstractHR{<:Complex}, v::Union{Tuple{<:Integer,<:Number},Number}, i, j, iRorR) = _HR_setindex!(hr, Val(true), v, i, j, iRorR)
"""
	numorb(obj) -> Int

Return the number of orbitals (or Wannier centres) represented by `obj`.
"""
@inline numorb(hr::AbstractHR) = hr.n
@inline Base.length(hr::AbstractHR) = length(hr.vals)
# @inline Base.size(hr::HR) = (hr.n, hr.n)
"""
	index(hr::AbstractHR, i::Integer) -> Tuple{Int, Int, Int}

Return the triple `(i, j, iR)` of the `i`-th stored term, where `iR` is the position of
the lattice vector in `hr.R`.
"""
@inline index(hr::AbstractHR, i::Integer) = (hr.is[i], hr.js[i], hr.iRs[i])
"""
	indices(hr::AbstractHR) -> iterator of (i, j, iR)

Iterate over the triples `(i, j, iR)` of the stored terms.
"""
@inline indices(hr::AbstractHR) = zip(hr.is, hr.js, hr.iRs)
@inline Base.enumerate(hr::AbstractHR) = zip(hr.is, hr.js, hr.iRs, hr.vals)
@inline Base.real(hr::AbstractHR{T}) where {T<:Real} = hr
@inline Base.complex(hr::AbstractHR{T}) where {T<:Complex} = hr
@inline Base.empty(::HT, n::Integer) where {HT<:AbstractHR} = empty(HT, n)
"""
	prune(hr::AbstractHR, cutoff::Real) -> AbstractHR

Return a copy of `hr` that keeps only the terms with `abs(v) ≥ cutoff` (for complex
values, `abs2(v) ≥ cutoff^2`).
"""
function prune(hr::AbstractHR{T}, cutoff::Real) where {T<:Real}
    @assert cutoff ≥ 0 "cutoff should be positive!"
    f(i, j, iR, v) = abs(v) ≥ cutoff
    return filter(f, hr)
end
function prune(hr::AbstractHR{T}, cutoff::Real) where {T<:Complex}
    @assert cutoff ≥ 0 "cutoff should be positive!"
    cutoff2 = cutoff^2
    f(i, j, iR, v) = abs2(v) ≥ cutoff2
    return filter(f, hr)
end
function _HR_sort_R(iRs, R, ndegen)
    perm = sortperm(R)
    R_new = R[perm]
    ndegen_new = ndegen[perm]
    perm_inv = invperm(perm)
    iRs_new = similar(iRs)
    @inbounds for (i, iR) in enumerate(iRs)
        iRs_new[i] = perm_inv[iR]
    end
    return iRs_new, R_new, ndegen_new
end
function _HR_sort_idx(is, js, iRs, vals)
    perm = sortperm(1:length(is); by=ih->(iRs[ih], js[ih], is[ih]))
    is_new = is[perm]
    js_new = js[perm]
    iRs_new = iRs[perm]
    vals_new = vals[perm]
    return is_new, js_new, iRs_new, vals_new
end
function _HR_sort_R!(iRs, R, ndegen)
    perm = sortperm(R)
    permute!(R, perm)
    permute!(ndegen, perm)
    perm_inv = invperm(perm)
    @inbounds for (i, iR) in enumerate(iRs)
        iRs[i] = perm_inv[iR]
    end
    return iRs, R, ndegen
end
function _HR_sort_idx!(is, js, iRs, vals)
    perm = sortperm(1:length(is); by=ih->(iRs[ih], js[ih], is[ih]))
    permute!(is, perm)
    permute!(js, perm)
    permute!(iRs, perm)
    permute!(vals, perm)
    return is, js, iRs, vals
end
"""
    HR(path, vals; storage = :full, kwargs...) -> HRFull or HRHermitian

Construct an `HR` object holding the data `v_{ijR}` given by the `(5, N)` matrix `path`
and the `N` values `vals`. Each column of `path` describes one term:
rows 1–3 are the lattice vector `R`, row 4 is the orbital index `i` and row 5 is the
orbital index `j`.

The keyword `storage` selects the representation: `:full` (default, see [`HRFull`](@ref))
stores every term independently, while `:hermi` (see [`HRHermitian`](@ref)) stores only
the independent terms of a Hermitian Hamiltonian, using `v_{ijR} = conj(v_{ji,-R})`.

# Example
```julia
julia> path = [ 1  0 -1;
                0  1  0;
                0  0  0;
                1  2  1;
                1  1  2 ];   # 5×3 matrix
julia> value = [1.0, 2.0, 3.0];
julia> hr = HR(path, value)
HR with 3 terms and 2 orbitals.
```
"""
HR(args...; storage::Symbol=:full, kwargs...) = _hr_from_val(Val(storage), args...; kwargs...)
_hr_from_val(::Val{:full}, args...; kwargs...) = HRFull(args...; kwargs...)
_hr_from_val(::Val{:hermi}, args...; kwargs...) = HRHermitian(args...; kwargs...)
_hr_from_val(::Val{S}, args...; kwargs...) where {S} =
    throw(ArgumentError("unknown storage: $S (expected :full or :hermi)"))

function _HR_filter_atol(is, js, iRsorRs, vals::AbstractVector{<:Complex}, atol)
    atol2 = atol^2
    fatol(v) = abs2(v) ≥ atol2
    mask = fatol.(vals)
    is = is[mask]
    js = js[mask]
    iRsorRs = iRsorRs[mask]
    vals = vals[mask]
    return is, js, iRsorRs, vals
end
function _HR_filter_atol(is, js, iRsorRs, vals::AbstractVector{<:Real}, atol)
    fatol(v) = abs(v) ≥ atol
    mask = fatol.(vals)
    is = is[mask]
    js = js[mask]
    iRsorRs = iRsorRs[mask]
    vals = vals[mask]
    return is, js, iRsorRs, vals
end

HRFull(hr::HRFull) = hr
function HRFull(hr::HRHermitian{T})::HRFull{T} where {T}
    nv = hr.nv
    double_nv = nv*2
    is = Vector{Int}(undef, 0)
    js = Vector{Int}(undef, 0)
    iRs = Vector{Int}(undef, 0)
    vals = Vector{Int}(undef, 0)
    sizehint!(is, double_nv)
    sizehint!(js, double_nv)
    sizehint!(iRs, double_nv)
    sizehint!(vals, double_nv)
    for ih in Base.OneTo(nv)
        ih2 = 2*ih
        ih1 = ih2 - 1
        i, j, iR, val = hr.is[ih], hr.js[ih], hr.iRs[ih], hr.vals[ih]
        push!(is, i)
        push!(js, j)
        push!(iRs, iR)
        push!(vals, val)
        if _HRHermitian_ispositive(i, j, iR, hr.R)
            push!(is, j)
            push!(js, i)
            push!(iRs, hr.negR[iR])
            push!(vals, conj(val))
        end
    end
    return HRFull{T}(hr.n, copy(hr.R), copy(hr.ndegen), is, js, iRs, vals, copy(hr.R2idx))
end
HRHermitian(hr::HRHermitian) = hr
function HRHermitian(hr::HRFull{T})::HRHermitian{T} where {T}
    negR = collect(get(hr.R2idx, -RR, 0) for RR in hr.R)
    any(iszero, negR) && throw(ArgumentError("hr may not be hermitian"))
    nv = hr.nv
    half_nv = floor(Int, nv/2) + numorb(hr) + 5
    unique_path2value = Dict{NTuple{3,Int},Tuple{Int,T}}()
    sizehint!(unique_path2value, half_nv)
    for ih in Base.OneTo(nv)
        i, j, iR, val = hr.is[ih], hr.js[ih], hr.iRs[ih], hr.vals[ih]
        if _HRHermitian_keep(i, j, iR, hr.R)
            path = (i, j, iR)
            (num, vsum) = get(unique_path2value, path, (0, zero(T)))
            unique_path2value[path] = (num+1, vsum + val)
        else
            path = (j, i, negR[iR])
            (num, vsum) = get(unique_path2value, path, (0, zero(T)))
            unique_path2value[path] = (num+1, vsum + conj(val))
        end
    end
    n = length(unique_path2value)
    is = Vector{Int}(undef, n)
    js = Vector{Int}(undef, n)
    iRs = Vector{Int}(undef, n)
    vals = Vector{T}(undef, n)
    i = 0
    for (key, value) in unique_path2value
        i += 1
        is[i] = key[1]
        js[i] = key[2]
        iRs[i] = key[3]
        vals[i] = value[2]
    end
    return HRHermitian{T}(hr.n, copy(hr.R), copy(hr.ndegen), is, js, iRs, vals, copy(hr.R2idx), negR)
end

# """
# 	reindexHR(hr::HR{T}, newindices::AbstractVector{<:Integer})::HR{T}

# This function will change the indices and path of hr.
# """
# function reindexHR(hr::HR{T}, newindices::AbstractVector{<:Integer})::HR{T} where {T}

#     if numorb(hr) ≠ length(newindices)
#         error("Wrong length of newindices to hr.")
#     end

#     neworbpath = Matrix{Int}(undef, 2, length(hr))
#     @views for i in eachindex(hr.indices)
#         I = hr.path[4:5, :] .== hr.indices[i]
#         neworbpath[I] .= newindices[i]
#     end

#     newpath = [hr.path[1:3, :]; neworbpath]

#     return HR{T}(newindices, newpath, copy(hr.value))
# end
