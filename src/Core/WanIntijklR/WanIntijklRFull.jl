struct WanIntijklRFull{T} <: WanIntijklR{T}
	nw::Int
	pR::Vector{ReducedCoordinates{Int}}
	pair::Vector{NTuple{3, Int}}
	R::Vector{ReducedCoordinates{Int}}
	ips::Vector{Int}
	jps::Vector{Int}
	iRs::Vector{Int}
	vals::Vector{T}
	pair2idx::Dict{NTuple{3, Int}, Int}
	pR2idx::Dict{ReducedCoordinates{Int}, Int}
	R2idx::Dict{ReducedCoordinates{Int}, Int}
end
function WanIntijklR_pairindex!(int::WanIntijklRFull, pair::WanIntPairIdxType)::Int
	ip = get(int.pair2idx, pair, 0)
	iszero(ip) || return ip
	npR = int.npR
	_isvalid_pair(pair, int.nw, npR) || error("pair $pair is invalid for nw=$(int.nw) and npR=$(npR)!")
	push!(int.pair, pair)
	ip = int.np
	int.pair2idx[pair] = ip
	return ip
end
function WanIntijklR_pairindex!(int::WanIntijklRFull, pair::WanIntPairVecType)::Int
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
		ip = WanIntijklR_pairindex!(int, pair)
	end
	return ip
end
function WanIntijklR_getindex(int::WanIntijklRFull{T}, ip1::Integer, ip2::Integer, iR::Integer) where {T}
	idx = findfirst(==((ip1, ip2, iR)), indices(int))
	return isnothing(idx) ? zero(T) : int.vals[idx]
end
function WanIntijklR_setindex!(int::WanIntijklRFull{T}, v::Number, ip1::Integer, ip2::Integer, iR::Integer) where {T}
	idx = findfirst(==((ip1, ip2, iR)), indices(int))
	isnothing(idx) || return int.vals[idx] = v
	push!(int.ips, ip1)
	push!(int.jps, ip2)
	push!(int.iRs, iR)
	push!(int.vals, v)
	return v
end
function Base.filter(f, int::WanIntijklRFull{T}) where {T}
	mask = collect(f(ip1, ip2, iR, v) for (ip1, ip2, iR, v) in enumerate(int))
	ips = int.ips[mask]
	jps = int.jps[mask]
	iRs = int.iRs[mask]
	vals = int.vals[mask]
	return WanIntijklRFull(ips, jps, iRs, vals, int.nw, int.pair, int.pR, int.R;
		inbound = true, dosort = false, dofilter = true)
end
@inline Base.zero(::Type{WanIntijklRFull{T}}, nw::Integer) where {T} =
	WanIntijklRFull{T}(nw,
		ReducedCoordinates{Int}[], NTuple{3, Int}[], ReducedCoordinates{Int}[],
		Int[], Int[], Int[], T[],
		Dict{NTuple{3, Int}, Int}(),
		Dict{ReducedCoordinates{Int}, Int}(),
		Dict{ReducedCoordinates{Int}, Int}())
@inline Base.real(int::WanIntijklRFull{T}) where {T <: Complex} =
	WanIntijklRFull{real(T)}(int.nw, copy(int.pR), copy(int.pair), copy(int.R),
		copy(int.ips), copy(int.jps), copy(int.iRs), real(int.vals),
		copy(int.pair2idx), copy(int.pR2idx), copy(int.R2idx))
function Base.show(io::IO, int::WanIntijklRFull)
	print(io, "WanIntijklR:")
	print(io, " nw: $(int.nw),")
	print(io, " np: $(int.np),")
	print(io, " nR: $(int.nR),")
	print(io, " nv: $(int.nv).")
end
"""
	_WanIntijklRFull_update!(ips, jps, iRs, pair, pR, R)

Detect the still existing pairs and R, and update ips, jps, iRs to the new indices. Also update pair and pR accordingly.
"""
function _WanIntijklRFull_update!(ips, jps, iRs, pair, pR, R)
	# pair
	pair_exist = sort(unique(vcat(ips, jps)))
	new_pair = pair[pair_exist]
	@inbounds if length(new_pair) < length(pair)
		pair_map = Dict(ip=>i for (i, ip) in enumerate(pair_exist))
		pair_map = [get(pair_map, ip, 0) for ip in eachindex(pair)]
		@simd for i in eachindex(ips) # @simd doesn't support enumerate.
			ips[i] = pair_map[ips[i]]
		end
		@simd for i in eachindex(jps)
			jps[i] = pair_map[jps[i]]
		end
	end
	# pR
	pR_exist = sort(unique(p[3] for p in new_pair))
	new_pR = pR[pR_exist]
	@inbounds if length(new_pR) < length(pR)
		pR_map = Dict(ipR=>i for (i, ipR) in enumerate(pR_exist))
		pR_map = [get(pR_map, ipR, 0) for ipR in eachindex(pR)]
		for (i, p) in enumerate(new_pair)
			new_pair[i] = (p[1], p[2], pR_map[p[3]])
		end
	end
	# R
	R_exist = sort(unique(iRs))
	new_R = R[R_exist]
	@inbounds if length(new_R) < length(R)
		R_map = Dict(iR=>i for (i, iR) in enumerate(R_exist))
		R_map = [get(R_map, iR, 0) for iR in eachindex(R)]
		@simd for i in eachindex(iRs)
			iRs[i] = R_map[iRs[i]]
		end
	end
	return ips, jps, iRs, new_pair, new_pR, new_R
end
function WanIntijklRFull(
	ips::AbstractVector{<:Integer},
	jps::AbstractVector{<:Integer},
	iRs::AbstractVector{<:Integer},
	vals::AbstractVector{T},
	nw::Integer,
	pair::AbstractVector{<:WanIntPairIdxType},
	pR::AbstractVector{<:AbstractVector{<:Integer}},
	R::AbstractVector{<:AbstractVector{<:Integer}};
	inbound::Bool = false, dosort::Bool = true, dofilter::Bool = true,
) where {T <: Number}

	if !inbound
		npR = length(pR)
		@assert all(p->_isvalid_pair(p, nw, npR), pair) "Invalid pair detected: nw=$(nw), npR=$(npR)."
		np = length(pair)
		@assert all(ip->1 ≤ ip ≤ np, ips) "Invalid ips detected: ip ∈ [1, $np]."
		@assert all(ip->1 ≤ ip ≤ np, jps) "Invalid jps detected: ip ∈ [1, $np]."
		nR = length(R)
		@assert all(iR->1 ≤ iR ≤ nR, iRs) "Invalid ips detected: iR ∈ [1, $nR]."
		@assert length(ips) == length(jps) == length(iRs) == length(vals) "Mismatched length of data."
	end

	isempty(vals) && return zero(WanIntijklRFull{T}, nw)

	pair = [(Int(p[1]), Int(p[2]), Int(p[3])) for p in pair]
	pR = ReducedCoordinates{Int}.(collect(pR))
	R = ReducedCoordinates{Int}.(collect(R))

	if dofilter
		# update data
		ips, jps, iRs, pair, pR, R = _WanIntijklRFull_update!(ips, jps, iRs, pair, pR, R)
	end

	if dosort
		# sort R by norm, and update iRs.
		pair, pR = _WanIntijklR_sortpR!(pair, pR)
		ips, jps, pair = _WanIntijklR_sortpair!(ips, jps, pair)
		iRs, R = _WanIntijklR_sortR!(iRs, R)
	end

	pair2idx = Dict(p=>i for (i, p) in enumerate(pair))
	pR2idx = Dict(RR=>i for (i, RR) in enumerate(pR))
	R2idx = Dict(RR=>i for (i, RR) in enumerate(R))

	return WanIntijklRFull{T}(nw, pR, pair, R,
		ips, jps, iRs, vals, pair2idx, pR2idx, R2idx)
end
function WanIntijklR(
	value::AbstractArray{T, 3},
	pair::AbstractVector{<:WanIntPairIdxType},
	pR::AbstractVector{<:AbstractVector{<:Integer}},
	R::AbstractVector{<:AbstractVector{<:Integer}};
	nw::Union{Integer, Nothing} = nothing, atol::Real = 1e-8,
	dosort::Bool = true, dofilter::Bool = true,
) where {T <: Number}
	@assert atol ≥ 0 "atol should be positive!"

	if isnothing(nw)
		nw = maximum(p->max(p[1], p[2]), pair)
	else
		nw = Int(nw)
	end

	npR = length(pR)
	nR = length(R)
	np = length(pair)
	@assert all(p->_isvalid_pair(p, nw, npR), pair) "Invalid pair detected: nw=$(nw), npR=$(npR)."

	size(value) == (np, np, nR) || error("Mismatched pair interaction!")

	ips = Int[]
	jps = Int[]
	iRs = Int[]
	vals = T[]

	estimated_len = nR * np * np  # 上界
	sizehint!(ips, estimated_len)
	sizehint!(jps, estimated_len)
	sizehint!(iRs, estimated_len)
	sizehint!(vals, estimated_len)

	fatol = _WanIntijklR_filter_atol(T, atol)
	@inbounds for iR in 1:nR, jp in 1:np, ip in 1:np
		v = value[ip, jp, iR]
		fatol(v) || continue
		push!(ips, ip)
		push!(jps, jp)
		push!(iRs, iR)
		push!(vals, v)
	end

	isempty(vals) && return zero(WanIntijklRFull{T}, nw)

	if T <: Complex && all(v -> abs(imag(v)) < atol, vals)
		vals = real(vals)
	end

	return WanIntijklRFull(ips, jps, iRs, vals, nw, pair, pR, R; inbound = true, dosort, dofilter)
end
function WanIntijklR(
	value::AbstractArray{T, 3},
	pair::AbstractVector{<:WanIntPairVecType},
	R::AbstractVector{<:AbstractVector{<:Integer}};
	kwargs...,
) where {T <: Number}
	pair_pR = collect(ReducedCoordinates{Int}(p[3]) for p in pair)
	pR = sort(unique(pair_pR); by = norm2)
	pR2idx = Dict{ReducedCoordinates{Int}, Int}(RR => ipR for (ipR, RR) in enumerate(pR))
	pair = map((ip, p)->(p[1], p[2], pR2idx[pair_pR[ip]]), enumerate(pair))
	return WanIntijklR(value, pair, pR, R; kwargs...)
end
function _WanIntijklR_filter_atol(::Type{<:Complex}, atol)
	atol2 = atol^2
	return v -> abs2(v) ≥ atol2
end
function _WanIntijklR_filter_atol(::Type{<:Real}, atol)
	return v -> abs(v) ≥ atol
end
