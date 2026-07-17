struct WanIntijklRHermitian{T} <: WanIntijklR{T}
	nw::Int
	pR::Vector{ReducedCoordinates{Int}}
	pair::Vector{NTuple{3, Int}} # all pairs
	R::Vector{ReducedCoordinates{Int}} # all R with conjR.
	# only contains half set of interactions.
	ips::Vector{Int}
	jps::Vector{Int}
	iRs::Vector{Int}
	vals::Vector{T}
	pair2idx::Dict{NTuple{3, Int}, Int}
	pR2idx::Dict{ReducedCoordinates{Int}, Int}
	R2idx::Dict{ReducedCoordinates{Int}, Int}
	# data for Hermitian
	uplo::Symbol # :U or :L
	uplopair::Set{Tuple{Symbol, Symbol}} # uolopairsym(uplo)
	conjpR::Vector{Int} # length(pR): pR[conjpR[ipR]] = -pR[ipR]
	pairtype::Vector{Symbol} # length(pair): :0, :p, :n; 
	pairpR::Vector{Int} # length(pair): pair[ip][3]; 
	conjpair::Vector{Int} # length(pair): pair[conjpair[ip]] = conj(pair[ip]), conj(p) = (p[2], p[1], -p[3])
	iconjRs::Vector{Int}
	# the set of R and pR define the set of conjR = R + pR2 - pR1.
	# Ridx::Set{Int}
	# conjRidx::Set{Int}
	# R′::Vector{ReducedCoordinates{Int}} # the set of pR - pR.
	# R′2idx::Dict{ReducedCoordinates{Int}, Int}
	# pR_minus_pR_to_R′::Ref{Matrix{Int}}
	# pair_minus_pair_to_R′::Ref{Matrix{Int}}
	# R_add_R′_to_R::Ref{Matrix{Int}} # (allR, R′), missing value is 0. If R + R′, then must have value, .
	# R_minus_R′_to_R::Ref{Matrix{Int}} # (allR, R′), missing value is 0. If conjR - R′, then must have value, .
end
include("./pre.jl")
include("./setindex!.jl")
function WanIntijklRHermitian_getindex(int::WanIntijklRHermitian{T}, ::Val{true}, ip1::Integer, ip2::Integer, iR::Integer) where {T}
	idx = findfirst(==((ip1, ip2, iR)), indices(int))
	return isnothing(idx) ? zero(T) : int.vals[idx]
end
function WanIntijklRHermitian_getindex(int::WanIntijklRHermitian{T}, ::Val{false}, ip1::Integer, ip2::Integer, iconjR::Integer) where {T}
	ip1 = int.conjpair[ip1]
	ip2 = int.conjpair[ip2]
	idx = findfirst(==((ip1, ip2, iconjR)), zip(int.ips, int.jps, int.iconjRs))
	return isnothing(idx) ? zero(T) : conj(int.vals[idx])
end
function WanIntijklR_getindex(int::WanIntijklRHermitian{T}, ip1::Integer, ip2::Integer, iR::Integer) where {T}
	p1t = int.pairtype[ip1]
	p2t = int.pairtype[ip2]
	return WanIntijklRHermitian_getindex(int, Val((p1t, p2t) ∈ int.uplopair), ip1, ip2, iR)
end
function Base.filter(f, int::WanIntijklRHermitian{T}) where {T}
	mask = collect(f(ip1, ip2, iR, v) for (ip1, ip2, iR, v) in enumerate(int))
	ips = int.ips[mask]
	jps = int.jps[mask]
	iRs = int.iRs[mask]
	vals = int.vals[mask]
	return WanIntijklRHermitian(ips, jps, iRs, vals,
		int.nw, int.pair, int.pR, int.R, int.uplo;
		ishermitian = true, inbound = true, dofilter = true)
end
@inline function Base.zero(::Type{WanIntijklRHermitian{T}}, nw::Integer, uplo::Symbol = :U) where {T}
	return WanIntijklRHermitian{T}(nw,
		ReducedCoordinates{Int}[], NTuple{3, Int}[], ReducedCoordinates{Int}[],
		Int[], Int[], Int[], T[],
		Dict{NTuple{3, Int}, Int}(),
		Dict{ReducedCoordinates{Int}, Int}(),
		Dict{ReducedCoordinates{Int}, Int}(),
		uplo, uplopairsym(Val(uplo)), Int[], Symbol[], Int[], Int[], Int[])
end
@inline Base.real(int::WanIntijklRHermitian{T}) where {T <: Complex} =
	WanIntijklRHermitian{real(T)}(int.nw,
		copy(int.pR), copy(int.pair), copy(int.R),
		copy(int.ips), copy(int.jps), copy(int.iRs), real(int.vals),
		copy(int.pair2idx), copy(int.pR2idx), copy(int.R2idx),
		int.uplo, copy(int.uplopair), copy(int.conjpR), copy(int.pairtype), copy(int.conjpair),
		copy(int.iconjRs))
function Base.show(io::IO, int::WanIntijklRHermitian)
	print(io, "WanIntijklR:")
	print(io, " nw: $(int.nw),")
	print(io, " np: $(int.np),")
	print(io, " nR: $(int.nR),")
	print(io, " nv: $(int.nv),")
	print(io, " uplo: $(int.uplo)")
end
function WanIntijklRHermitian(
	ips::AbstractVector{<:Integer},
	jps::AbstractVector{<:Integer},
	iRs::AbstractVector{<:Integer},
	vals::AbstractVector{T},
	nw::Integer,
	pair::AbstractVector{<:WanIntPairIdxType},
	pR::AbstractVector{<:AbstractVector{<:Integer}},
	R::AbstractVector{<:AbstractVector{<:Integer}},
	uplo::Symbol; ishermitian::Bool = false,
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

	pair = [(Int(p[1]), Int(p[2]), Int(p[3])) for p in pair]
	pR = ReducedCoordinates{Int}.(collect(pR))
	R = ReducedCoordinates{Int}.(collect(R))

	if dofilter
		ips, jps, iRs, pair, pR, R = _WanIntijklRFull_update!(ips, jps, iRs, pair, pR, R)
	end

	# preprocessing is light.
	(ips, jps, pair, pair2idx, conjpair, pR, pR2idx, conjpR) = WanIntijklRHermitian_pre_pair!(ips, jps, pair, pR)
	pairtype = typeofpair(pair, pR)
	pairpR = [p[3] for p in pair]
	uplopair = uplopairsym(Val(uplo))
	if !ishermitian
		(ips, jps, iRs, vals, R) = WanIntijklRHermitian_pre_value(ips, jps, iRs, vals, pR, R, conjpair, pairtype, pairpR, uplopair)
	end
	(iRs, iconjRs, R, R2idx) = WanIntijklRHermitian_pre_iconjR!(ips, jps, iRs, pair, pR, R, pairtype, uplopair)

	return WanIntijklRHermitian{T}(nw, pR, pair, R,
		ips, jps, iRs, vals, pair2idx, pR2idx, R2idx,
		uplo, uplopair, conjpR, pairtype, pairpR, conjpair,
		iconjRs)
end
function WanIntijklR(
	value::AbstractArray{T, 3},
	pair_idx::AbstractVector{<:NTuple{2, <:Integer}},
	pair::AbstractVector{<:WanIntPairIdxType},
	pR::AbstractVector{<:AbstractVector{<:Integer}},
	R::AbstractVector{<:AbstractVector{<:Integer}};
	nw::Union{Integer, Nothing} = nothing, atol::Real = 1e-8,
	uplo::Union{Symbol, Nothing} = nothing, ishermitian::Bool = true,
	dosort::Bool = true, dofilter::Bool = true,
) where {T <: Number}
	@assert atol ≥ 0 "atol should be positive!"

	# 数据只读pair_idx; ishermitian决定了后续WanIntijklRHermitian是否要对源数据进行厄密化处理。
	if isnothing(nw)
		nw = maximum(p->max(p[1], p[2]), pair)
	else
		nw = Int(nw)
	end

	npR = length(pR)
	nR = length(R)
	np = length(pair)
	@assert all(p->_isvalid_pair(p, nw, npR), pair) "Invalid pair detected: nw=$(nw), npR=$(npR)."

	if isnothing(uplo)
		pairtype = typeofpair(pair, pR)
		pp = pair_idx[1]
		pptype = (pairtype[pp[1]], pairtype[pp[2]])
		if pptype ∈ uplopairsym(Val(:U))
			uplo = :U
		else
			uplo = :L
		end
	elseif uplo ∉ [:U, :L]
		error("uplo should be :U or :L!")
	end

	size(value) == (np, np, nR) || error("Mismatched pair interaction!")

	ips = Int[]
	jps = Int[]
	iRs = Int[]
	vals = T[]

	estimated_len = ceil(Int, nR * np * np / 2)  # 上界
	sizehint!(ips, estimated_len)
	sizehint!(jps, estimated_len)
	sizehint!(iRs, estimated_len)
	sizehint!(vals, estimated_len)

	fatol = _WanIntijklR_filter_atol(T, atol)
	@inbounds for iR in 1:nR, (ip, jp) in pair_idx
		v = value[ip, jp, iR]
		fatol(v) || continue
		push!(ips, ip)
		push!(jps, jp)
		push!(iRs, iR)
		push!(vals, v)
	end

	isempty(vals) && return zero(WanIntijklRHermitian{T}, nw)

	if T <: Complex && all(v -> abs(imag(v)) < atol, vals)
		vals = real(vals)
	end

	return WanIntijklRHermitian(ips, jps, iRs, vals,
		nw, pair, pR, R, uplo;
		ishermitian, inbound = true, dosort, dofilter)
end
function WanIntijklR(
	value::AbstractArray{T, 3},
	pair_idx::AbstractVector{<:NTuple{2, <:Integer}},
	pair::AbstractVector{<:WanIntPairVecType},
	R::AbstractVector{<:AbstractVector{<:Integer}};
	kwargs...,
) where {T <: Number}
	pair_pR = collect(ReducedCoordinates{Int}(p[3]) for p in pair)
	pR = sort(unique(pair_pR); by = norm2)
	pR2idx = Dict{ReducedCoordinates{Int}, Int}(RR => ipR for (ipR, RR) in enumerate(pR))
	pair = map((ip, p)->(p[1], p[2], pR2idx[pair_pR[ip]]), enumerate(pair))
	return WanIntijklR(value, pair_idx, pair, pR, R; kwargs...)
end
