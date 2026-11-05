export PairIntHR, PairIntHRFull, PairIntHRHermitian
export numorb, index, indices, prune, filter_pair
"""
	PairIntHR{T<:Number}

Abstract supertype of pair interaction matrix elements ``V_{p_1p_2R}``, indexed by
orbital pairs `p` and lattice vectors `R`.
"""
abstract type PairIntHR{T<:Number} end
# PairIntHR
# 底层数据存储依赖于AbstractHR
# 支持访问和修改元素，但是该操作较为低效，仅用作检查用。
function Base.getproperty(int::PairIntHR, name::Symbol)
    if name === :npR
        return length(getfield(int, :pR))
    elseif name === :np
        return length(getfield(int, :pair))
    elseif name === :nR
        return getfield(int, :hr).nR
    elseif name === :nv
        return length(getfield(int, :hr))
    elseif name === :ips
        return getfield(getfield(int, :hr), :is)
    elseif name === :jps
        return getfield(getfield(int, :hr), :js)
    elseif name === :iRs
        return getfield(getfield(int, :hr), :iRs)
    elseif name === :vals
        return getfield(getfield(int, :hr), :vals)
    elseif name === :R
        return getfield(getfield(int, :hr), :R)
    elseif name === :ndegen
        return getfield(getfield(int, :hr), :ndegen)
    elseif name === :negR
        return getfield(getfield(int, :hr), :negR)
    else
        return getfield(int, name)
    end
end
function Base.show(io::IO, int::PairIntHR)
    print(io, "PairIntHR:")
    print(io, " nw: $(int.nw),")
    print(io, " np: $(int.np),")
    print(io, " nR: $(int.nR),")
    print(io, " nv: $(int.nv).")
end

include("./pair.jl")
include("./Full.jl")
include("./Hermitian.jl")
# include("./generate.jl")

function PairIntHR_pairindex(int::PairIntHR, ip::Integer)::Int
    np = int.np
    1 ≤ ip ≤ np || error("ip = $ip should be in [1, $np]")
    return ip
end
function PairIntHR_pairindex(int::PairIntHR, pair::WanIntPairIdxType)::Int
    return get(int.pair2idx, pair, 0)
end
function PairIntHR_pairindex(int::PairIntHR, pair::WanIntPairVecType)::Int
    pR = ReducedCoordinates{Int}(pair[3])
    ipR = get(int.pR2idx, pR, 0)
    iszero(ipR) && return 0
    pair = (pair[1], pair[2], ipR)
    return PairIntHR_pairindex(int, pair)
end
PairIntHR_pairindex!(int::PairIntHR, ip::Integer) = PairIntHR_pairindex(int, ip)
function PairIntHR_pairindex!(int::PairIntHR, pair::WanIntPairIdxType)::Int
    ip = get(int.pair2idx, pair, 0)
    iszero(ip) || return ip
    _isvalid_pair(pair, int.nw, int.npR) || error("pair $pair is invalid for nw=$(int.nw) and npR=$(int.npR)!")
    push!(int.pair, pair)
    ip = int.np
    int.pair2idx[pair] = ip
    return ip
end
function PairIntHR_pairindex!(int::PairIntHR, pair::WanIntPairVecType)::Int
    pR = ReducedCoordinates{Int}(pair[3])
    ipR = get(int.pR2idx, pR, 0)
    if iszero(ipR)
        _isvalid_pair(pair, int.nw) || error("pair $pair is invalid for nw=$(int.nw)!")
        push!(int.pR, pR)
        ipR = int.npR
        int.pR2idx[pR] = ipR
        pair = (pair[1], pair[2], ipR)
        push!(int.pair, pair)
        ip = int.np
        int.pair2idx[pair] = ip
    else
        pair = (pair[1], pair[2], ipR)
        ip = PairIntHR_pairindex!(int, pair)
    end
    return ip
end
@inline Base.getindex(int::PairIntHR, i::Integer) = int.hr[i]
function Base.getindex(int::PairIntHR{T}, ip1orp1, ip2orp2, iRorR) where {T}
    ip1 = PairIntHR_pairindex(int, ip1orp1)
    iszero(ip1) && return zero(T)
    ip2 = PairIntHR_pairindex(int, ip2orp2)
    iszero(ip2) && return zero(T)
    return _HR_getindex(int.hr, Val(false), ip1, ip2, iRorR)
end
@inline Base.setindex!(int::PairIntHR{T}, v::Real, i::Integer) where {T<:Real} = (int.hr[i] = v)
@inline Base.setindex!(int::PairIntHR{T}, v::Number, i::Integer) where {T<:Complex} = (int.hr[i] = v)
function Base.setindex!(int::PairIntHR{T}, v::Union{Tuple{<:Integer,<:Real},Real}, ip1orp1, ip2orp2, iRorR) where {T<:Real}
    ip1 = PairIntHR_pairindex!(int, ip1orp1)
    ip2 = PairIntHR_pairindex!(int, ip2orp2)
    return _HR_setindex!(int.hr, Val(false), v, ip1, ip2, iRorR)
end
function Base.setindex!(int::PairIntHR{T}, v::Union{Tuple{<:Integer,<:Number},Number}, ip1orp1, ip2orp2, iRorR) where {T<:Complex}
    ip1 = PairIntHR_pairindex!(int, ip1orp1)
    ip2 = PairIntHR_pairindex!(int, ip2orp2)
    return _HR_setindex!(int.hr, Val(false), v, ip1, ip2, iRorR)
end
@inline numorb(int::PairIntHR) = int.nw
@inline Base.length(int::PairIntHR) = length(int.hr)
@inline index(int::PairIntHR, i::Integer) = index(int.hr, i)
@inline indices(int::PairIntHR) = indices(int.hr)
@inline Base.enumerate(int::PairIntHR) = enumerate(int.hr)
@inline Base.real(int::PairIntHR{T}) where {T<:Real} = int
@inline Base.empty(::T, args...) where {T<:PairIntHR} = empty(T, args...)
function prune(int::PairIntHR{T}, cutoff::Real) where {T<:Real}
    @assert cutoff ≥ 0 "cutoff should be positive!"
    f(ip1, ip2, iR, v) = abs(v) ≥ cutoff
    return filter(f, int)
end
function prune(int::PairIntHR{T}, cutoff::Real) where {T<:Complex}
    @assert cutoff ≥ 0 "cutoff should be positive!"
    cutoff2 = cutoff^2
    f(ip1, ip2, iR, v) = abs2(v) ≥ cutoff2
    return filter(f, int)
end
"""
	filter_pair(int::PairIntHR, pair::AbstractVector{<:WanIntPairVecType}) -> PairIntHR

Return the sub-block of `int` that involves only the orbital pairs listed in `pair`,
given as `(iw1, iw2, R)` tuples. If no stored pair is selected, return an empty object of
the same type and `nw`.
"""
function filter_pair(int::PairIntHR, pair::AbstractVector{<:WanIntPairVecType})

    pair = Set{Tuple{Int,Int,ReducedCoordinates{Int}}}(pair)
    isempty(pair) && return zero(int, int.nw)

    hold_pair_idx = findall(int.pair) do (iw1, iw2, ipR)
        (iw1, iw2, int.pR[ipR]) ∈ pair
    end
    isempty(hold_pair_idx) && return zero(int, int.nw)

    f(ip1, ip2, iR, v) = (ip1 ∈ hold_pair_idx && ip2 ∈ hold_pair_idx)

    return filter(f, int)
end
function _PairIntHR_sortpR!(pair, pR)
    pRsortidx = sortperm(pR; by=norm2)
    pR = pR[pRsortidx]
    pRsortidx_inv = invperm(pRsortidx)
    @inbounds for (i, p) in enumerate(pair)
        pair[i] = (p[1], p[2], pRsortidx_inv[p[3]])
    end
    return pair, pR
end
function _PairIntHR_sortpair!(ips, jps, pair)
    by(p) = (p[3], p[1], p[2])
    pairsortidx = sortperm(pair; by)
    pair = pair[pairsortidx]
    pairsortidx_inv = invperm(pairsortidx)
    @inbounds for (i, ip) in enumerate(ips)
        ips[i] = pairsortidx_inv[ip]
    end
    @inbounds for (i, ip) in enumerate(jps)
        jps[i] = pairsortidx_inv[ip]
    end
    return ips, jps, pair
end
function _PairIntHR_sortR!(iRs, R)
    Rsortidx = sortperm(R; by=norm2)
    R = R[Rsortidx]
    Rsortidx_inv = invperm(Rsortidx)
    @inbounds for (i, iR) in enumerate(iRs)
        iRs[i] = Rsortidx_inv[iR]
    end
    return iRs, R
end
"""
	PairIntHR(args...; storage::Symbol=:full, kwargs...)

Construct a pair interaction container from the given data. The keyword `storage`
selects [`PairIntHRFull`](@ref) (`:full`) or [`PairIntHRHermitian`](@ref) (`:hermi`).
"""
PairIntHR(args...; storage::Symbol=:full, kwargs...) = _PairIntHR_from_val(Val(storage), args...; kwargs...)
_PairIntHR_from_val(::Val{:full}, args...; kwargs...) = PairIntHRFull(args...; kwargs...)
_PairIntHR_from_val(::Val{:hermi}, args...; kwargs...) = PairIntHRHermitian(args...; kwargs...)
_PairIntHR_from_val(::Val{S}, args...; kwargs...) where {S} =
    throw(ArgumentError("unknown storage: $S (expected :full or :hermi)"))
function _PairIntHR_filter_idx!(ips, jps, pair)
    # pair
    pair_exist = unique(vcat(ips, jps))
    sort!(pair_exist)
    new_pair = unique(view(pair, pair_exist))
    @inbounds if length(new_pair) < length(pair)
        new_pair_map = Dict{Tuple{Int,Int,ReducedCoordinates{Int}},Int}(p=>i for (i, p) in enumerate(new_pair))
        index_to_new = [get(new_pair_map, (p[1], p[2], pR[p[3]]), 0) for p in pair]
        for (i, ip) in enumerate(ips) # @simd doesn't support enumerate.
            ips[i] = index_to_new[ip]
        end
        for (i, ip) in enumerate(jps)
            jps[i] = index_to_new[ip]
        end
    else
        new_pair = pair
    end
    return ips, jps, new_pair
end
function _PairIntHR_filter_idx!(ips, jps, pair, pR)
    # pair
    pair_exist = unique(vcat(ips, jps))
    sort!(pair_exist)
    new_pair = map(pair_exist) do ip
        ip, jp, ipR = pair[ip]
        return (ip, jp, pR[ipR])
    end
    unique!(new_pair)
    @inbounds if length(new_pair) < length(pair)
        new_pair_map = Dict{Tuple{Int,Int,ReducedCoordinates{Int}},Int}(p=>i for (i, p) in enumerate(new_pair))
        index_to_new = [get(new_pair_map, (p[1], p[2], pR[p[3]]), 0) for p in pair]
        for (i, ip) in enumerate(ips) # @simd doesn't support enumerate.
            ips[i] = index_to_new[ip]
        end
        for (i, ip) in enumerate(jps)
            jps[i] = index_to_new[ip]
        end
        pR2idx = Dict{ReducedCoordinates{Int},Int}(RR=>i for (i, RR) in enumerate(pR))
        new_pair = [(p[1], p[2], pR2idx[p[3]]) for p in new_pair]
    else
        new_pair = pair
    end
    # pR
    pR_exist = unique(p[3] for p in new_pair)
    sort!(pR_exist)
    new_pR = unique(view(pR, pR_exist))
    @inbounds if length(new_pR) < length(pR)
        new_pR_map = Dict{ReducedCoordinates{Int},Int}(RR=>i for (i, RR) in enumerate(new_pR))
        index_to_new = [get(new_pR_map, RR, 0) for RR in pR]
        for (i, p) in enumerate(new_pair)
            new_pair[i] = (p[1], p[2], index_to_new[p[3]])
        end
    else
        new_pR = pR
    end
    return ips, jps, new_pair, new_pR
end

PairIntHRFull(int::PairIntHRFull) = int
PairIntHRFull(int::PairIntHRHermitian{T}) where {T} =
    PairIntHRFull{T}(int.nw, int.pR, int.pair, HRFull(int.hr), int.pR2idx, int.pair2idx)
PairIntHRHermitian(int::PairIntHRHermitian) = int
PairIntHRHermitian(int::PairIntHRFull{T}) where {T} =
    PairIntHRHermitian{T}(int.nw, int.pR, int.pair, HRHermitian(int.hr), int.pR2idx, int.pair2idx)
