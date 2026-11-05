"""
	PairIntHRHermitian{T}

Hermitian storage of pair interaction matrix elements: only the independent terms under
`V_{p_1p_2R} = conj(V_{p_2p_1,-R})` are kept.
"""
struct PairIntHRHermitian{T} <: PairIntHR{T}
    nw::Int
    pR::Vector{ReducedCoordinates{Int}}
    pair::Vector{NTuple{3,Int}}
    hr::HRHermitian{T}
    pR2idx::Dict{ReducedCoordinates{Int},Int}
    pair2idx::Dict{NTuple{3,Int},Int}
end
@inline Base.real(int::PairIntHRHermitian{T}) where {T<:Complex} =
    PairIntHRHermitian{real(T)}(int.nw, copy(int.pR), copy(int.pair), real(int.hr), copy(int.pR2idx), copy(int.pair2idx))
@inline Base.empty(::Type{PairIntHRHermitian{T}}, nw::Integer) where {T} =
    PairIntHRHermitian{T}(nw, ReducedCoordinates{Int}[], NTuple{3,Int}[], empty(HRHermitian{T}, typemax(Int)),
        Dict{ReducedCoordinates{Int},Int}(), Dict{NTuple{3,Int},Int}())
@inline Base.convert(::Type{PairIntHRHermitian{T₁}}, int::PairIntHRHermitian{T₂}) where {T₁,T₂} =
    PairIntHRHermitian{T₁}(int.nw, copy(int.pR), copy(int.pair), convert(HRHermitian{T₁}, int.hr),
        copy(int.pR2idx), copy(int.pair2idx))
function Base.filter(f, int::PairIntHRHermitian{T}) where {T}
    mask = collect(f(ip1, ip2, iR, v) for (ip1, ip2, iR, v) in enumerate(int))
    ips = int.ips[mask]
    jps = int.jps[mask]
    iRs = int.iRs[mask]
    vals = int.vals[mask]
    return PairIntHRHermitian(ips, jps, iRs, vals, int.nw, int.pair, int.pR, int.R; ndegen=int.ndegen, inbound=true, dofilter=true)
end
# function Base.sort(int::PairIntHRHermitian{T}, lattice, orblocat) where {T}
# end
function PairIntHRHermitian(
    ips::AbstractVector{<:Integer},
    jps::AbstractVector{<:Integer},
    iRs::AbstractVector{<:Integer},
    vals::AbstractVector{T},
    nw::Integer,
    pair::AbstractVector{<:WanIntPairVecType},
    R::AbstractVector{<:AbstractVector{<:Integer}};
    ndegen::AbstractVector{<:Integer}=ones(Int, length(R)),
    inbound::Bool=false, dofilter::Bool=true, have_negR::Bool=false, docheck::Bool=true,
) where {T<:Number}

    if !inbound
        @assert all(p->_isvalid_pair(p, nw), pair) "Invalid pair detected: nw=$(nw)."
        np = length(pair)
        @assert all(ip->1 ≤ ip ≤ np, ips) "Invalid ips detected: ip ∈ [1, $np]."
        @assert all(ip->1 ≤ ip ≤ np, jps) "Invalid jps detected: ip ∈ [1, $np]."
        nR = length(R)
        @assert all(iR->1 ≤ iR ≤ nR, iRs) "Invalid ips detected: iR ∈ [1, $nR]."
        @assert length(ndegen) == nR "Mismatched length of ndegen."
        @assert all(>(0), ndegen) "ndegen should be positive."
        @assert length(ips) == length(jps) == length(iRs) == length(vals) "Mismatched length of data."
    end
    isempty(vals) && return empty(PairIntHRHermitian{T}, nw)

    pair = [(Int(p[1]), Int(p[2]), ReducedCoordinates{Int}(p[3])) for p in pair]
    R = collect(ReducedCoordinates{Int}, R)
    ndegen = collect(Int, ndegen)

    if dofilter
        ips, jps, pair = _PairIntHR_filter_idx!(ips, jps, pair)
    end

    pR = unique(p[3] for p in pair)
    pR2idx = Dict(RR=>i for (i, RR) in enumerate(pR))
    pair = [(p[1], p[2], pR2idx[p[3]]) for p in pair]
    pair2idx = Dict(p=>i for (i, p) in enumerate(pair))

    hr = HRHermitian(ips, jps, iRs, vals, typemax(Int), R; ndegen, inbound=true, dofilter, have_negR, docheck)
    return PairIntHRHermitian{T}(nw, pR, pair, hr, pR2idx, pair2idx)
end
function PairIntHRHermitian(
    ips::AbstractVector{<:Integer},
    jps::AbstractVector{<:Integer},
    iRs::AbstractVector{<:Integer},
    vals::AbstractVector{T},
    nw::Integer,
    pair::AbstractVector{<:WanIntPairIdxType},
    pR::AbstractVector{<:AbstractVector{<:Integer}},
    R::AbstractVector{<:AbstractVector{<:Integer}};
    ndegen::AbstractVector{<:Integer}=ones(Int, length(R)),
    inbound::Bool=false, dofilter::Bool=true, have_negR::Bool=false, docheck::Bool=true,
) where {T<:Number}

    if !inbound
        npR = length(pR)
        @assert all(p->_isvalid_pair(p, nw, npR), pair) "Invalid pair detected: nw=$(nw), npR=$(npR)."
        np = length(pair)
        @assert all(ip->1 ≤ ip ≤ np, ips) "Invalid ips detected: ip ∈ [1, $np]."
        @assert all(ip->1 ≤ ip ≤ np, jps) "Invalid jps detected: ip ∈ [1, $np]."
        nR = length(R)
        @assert all(iR->1 ≤ iR ≤ nR, iRs) "Invalid ips detected: iR ∈ [1, $nR]."
        @assert length(ndegen) == nR "Mismatched length of ndegen."
        @assert all(>(0), ndegen) "ndegen should be positive."
        @assert length(ips) == length(jps) == length(iRs) == length(vals) "Mismatched length of data."
    end
    isempty(vals) && return empty(PairIntHRHermitian{T}, nw)

    pair = [(Int(p[1]), Int(p[2]), Int(p[3])) for p in pair]
    pR = collect(ReducedCoordinates{Int}, pR)
    R = collect(ReducedCoordinates{Int}, R)
    ndegen = collect(Int, ndegen)

    if dofilter
        ips, jps, pair, pR = _PairIntHR_filter_idx!(ips, jps, pair, pR)
    end

    pair2idx = Dict(p=>i for (i, p) in enumerate(pair))
    pR2idx = Dict(RR=>i for (i, RR) in enumerate(pR))

    hr = HRHermitian(ips, jps, iRs, vals, typemax(Int), R; ndegen, inbound=true, dofilter, have_negR, docheck)
    return PairIntHRHermitian{T}(nw, pR, pair, hr, pR2idx, pair2idx)
end
