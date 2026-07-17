export WanIntijklRFull, WanIntijklRHermitian,
	numorb, index, indices, prune, filter_pair

# struct <: WanIntijklR{T}
# 	nw::Int
# 	pR::Vector{ReducedCoordinates{Int}}
# 	pair::Vector{NTuple{3, Int}} # (iw1, iw2, ipR)
# 	R::Vector{ReducedCoordinates{Int}}
# 	ips::Vector{Int}
# 	jps::Vector{Int}
# 	iRs::Vector{Int}
# 	vals::Vector{T}
# 	# dict used to getindex and setindex!
# 	pair2idx::Dict{NTuple{3, Int}, Int}
# 	pR2idx::Dict{ReducedCoordinates{Int}, Int}
# 	R2idx::Dict{ReducedCoordinates{Int}, Int}
# end
# 支持访问和修改元素，但是该操作极为低效，仅用作检查用；
# 推荐访问方式为遍历。

function Base.getproperty(int::WanIntijklR, name::Symbol)
	if name === :npR
		return length(int.pR)
	elseif name === :np
		return length(int.pair)
	elseif name === :nR
		return length(int.R)
	elseif name === :nv
		return length(int.vals)
	else
		return getfield(int, name)
	end
end

include("./pair.jl")
include("./WanIntijklRFull.jl")
include("./WanIntijklRHermitian/WanIntijklRHermitian.jl")

function WanIntijklR_pairindex(int::WanIntijklR, ip::Integer)::Int
	np = int.np
	1 ≤ ip ≤ np || error("ip = $ip should be in [1, $np]")
	return ip
end
function WanIntijklR_pairindex(int::WanIntijklR, pair::WanIntPairIdxType)::Int
	return get(int.pair2idx, pair, 0)
end
function WanIntijklR_pairindex(int::WanIntijklR, pair::WanIntPairVecType)::Int
	pR = ReducedCoordinates{Int}(pair[3])
	ipR = get(int.pR2idx, pR, 0)
	iszero(ipR) && return 0
	pair = (pair[1], pair[2], ipR)
	return WanIntijklR_pairindex(int, pair)
end
function WanIntijklR_Rindex(int::WanIntijklR, iR::Integer)::Int
	nR = int.nR
	1 ≤ iR ≤ nR || error("iR = $iR should be in [1, $nR]")
	return iR
end
function WanIntijklR_Rindex(int::WanIntijklR, R::AbstractVector{<:Integer})::Int
	R = ReducedCoordinates{Int}(R)
	return get(int.R2idx, R, 0)
end
WanIntijklR_pairindex!(int::WanIntijklR, ip::Integer) = WanIntijklR_pairindex(int, ip)
WanIntijklR_Rindex!(int::WanIntijklR, iR::Integer) = WanIntijklR_Rindex(int, iR)
function WanIntijklR_Rindex!(int::WanIntijklR, R::AbstractVector{<:Integer})::Int
	R = ReducedCoordinates{Int}(R)
	iR = get(int.R2idx, R, 0)
	if iszero(iR)
		push!(int.R, R)
		iR = int.nR
		int.R2idx[R] = iR
	end
	return iR
end

@inline Base.getindex(int::WanIntijklR, i::Integer) = int.vals[i]
function Base.getindex(int::WanIntijklR{T}, ip1orp1, ip2orp2, iRorR) where {T}
	ip1 = WanIntijklR_pairindex(int, ip1orp1)
	iszero(ip1) && return zero(T)
	ip2 = WanIntijklR_pairindex(int, ip2orp2)
	iszero(ip2) && return zero(T)
	iR = WanIntijklR_Rindex(int, iRorR)
	iszero(iR) && return zero(T)
	return WanIntijklR_getindex(int, ip1, ip2, iR)
end
@inline Base.setindex!(int::WanIntijklR{<:Real}, v::Real, i::Integer) = (int.vals[i] = v)
@inline Base.setindex!(int::WanIntijklR{<:Complex}, v::Number, i::Integer) = (int.vals[i] = v)
function Base.setindex!(int::WanIntijklR{<:Real}, v::Real, ip1orp1, ip2orp2, iRorR)
	ip1 = WanIntijklR_pairindex!(int, ip1orp1)
	ip2 = WanIntijklR_pairindex!(int, ip2orp2)
	iR = WanIntijklR_Rindex!(int, iRorR)
	return WanIntijklR_setindex!(int, v, ip1, ip2, iR)
end
function Base.setindex!(int::WanIntijklR{<:Complex}, v::Number, ip1orp1, ip2orp2, iRorR)
	ip1 = WanIntijklR_pairindex!(int, ip1orp1)
	ip2 = WanIntijklR_pairindex!(int, ip2orp2)
	iR = WanIntijklR_Rindex!(int, iRorR)
	return WanIntijklR_setindex!(int, v, ip1, ip2, iR)
end

@inline numorb(int::WanIntijklR) = int.nw
@inline Base.length(int::WanIntijklR) = length(int.vals)
@inline index(int::WanIntijklR, i::Integer) = (int.ips[i], int.jps[i], int.iRs[i])
@inline indices(int::WanIntijklR) = zip(int.ips, int.jps, int.iRs)
@inline Base.enumerate(int::WanIntijklR) = zip(int.ips, int.jps, int.iRs, int.vals)
@inline Base.real(int::WanIntijklR{T}) where {T <: Real} = int
@inline Base.zero(::T, args...) where {T <: WanIntijklR} = zero(T, args...)
function prune(int::WanIntijklR{T}, cutoff::Real) where {T <: Real}
	@assert cutoff ≥ 0 "cutoff should be positive!"
	f(ip1, ip2, iR, v) = abs(v) ≥ cutoff
	return filter(f, int)
end
function prune(int::WanIntijklR{T}, cutoff::Real) where {T <: Complex}
	@assert cutoff ≥ 0 "cutoff should be positive!"
	cutoff2 = cutoff^2
	f(ip1, ip2, iR, v) = abs2(v) ≥ cutoff2
	return filter(f, int)
end
function filter_pair(int::WanIntijklR, pair::AbstractVector)

	pair = Set{NTuple{5, Int}}(pair)
	isempty(pair) && return zero(int, int.nw)

	hold_pair_idx = findall(int.pair) do (iw1, iw2, ipR)
		pR = int.pR[ipR]
		(iw1, iw2, pR[1], pR[2], pR[3]) ∈ pair
	end
	isempty(hold_pair_idx) && return zero(int, int.nw)

	f(ip1, ip2, iR, v) = (ip1 ∈ hold_pair_idx && ip2 ∈ hold_pair_idx)

	return filter(f, int)
end
#TODO sort
function hermitianpart(int::WanIntijklRFull{T}) where {T}
	return WanIntijklRHermitian(
		int.ips, int.jps, int.iRs, int.vals,
		int.nw, int.pair, int.pR, int.R,
		:U; ishermitian = false, inbound = true, dofilter = true)
end

function _WanIntijklR_sortpR!(pair, pR)
	pRsortidx = sortperm(pR; by = norm2)
	pR = pR[pRsortidx]
	pRsortidx_inv = invperm(pRsortidx)
	@inbounds for (i, p) in enumerate(pair)
		pair[i] = (p[1], p[2], pRsortidx_inv[p[3]])
	end
	return pair, pR
end
function _WanIntijklR_sortpair!(ips, jps, pair)
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
function _WanIntijklR_sortR!(iRs, R)
	Rsortidx = sortperm(R; by = norm2)
	R = R[Rsortidx]
	Rsortidx_inv = invperm(Rsortidx)
	@inbounds for (i, iR) in enumerate(iRs)
		iRs[i] = Rsortidx_inv[iR]
	end
	return iRs, R
end
