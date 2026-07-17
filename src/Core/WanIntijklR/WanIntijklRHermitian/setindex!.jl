function WanIntijklRHermitian_pRindex!(int::WanIntijklRHermitian, pR::ReducedCoordinates{Int})::Int
	ipR = get(int.pR2idx, pR, 0)
	iszero(ipR) || return ipR
	push!(int.pR, pR)
	push!(int.pR, -pR)
	npR = int.npR
	ipR = npR - 1
	int.pR2idx[pR] = ipR
	int.pR2idx[-pR] = npR
	push!(int.conjpR, npR)
	push!(int.conjpR, ipR)
	return ipR
end
function WanIntijklR_pairindex!(int::WanIntijklRHermitian, pair::WanIntPairIdxType)::Int
	ip = get(int.pair2idx, pair, 0)
	iszero(ip) || return ip
	npR = int.npR
	_isvalid_pair(pair, int.nw, npR) || error("pair $pair is invalid for nw=$(int.nw) and npR=$(npR)!")
	push!(int.pair, pair)
	pt = typeofpair(pair, int.pR)
	push!(int.pairtype, pt)
	push!(int.pairpR, pair[3])
	ip = int.np
	int.pair2idx[pair] = ip
	if pt == Symbol(0)
		push!(int.conjpair, ip)
	else
		conjpair = (pair[2], pair[1], int.conjpR[pair[3]])
		push!(int.pair, conjpair)
		push!(int.pairtype, typeofpair(conjpair, int.pR))
		push!(int.pairpR, conjpair[3])
		np = int.np
		int.pair2idx[conjpair] = np
		push!(int.conjpair, np)
		push!(int.conjpair, ip)
	end
	return ip
end
function WanIntijklR_pairindex!(int::WanIntijklRHermitian, pair::WanIntPairVecType)::Int
	pR = ReducedCoordinates{Int}(pR)
	ipR = WanIntijklRHermitian_pRindex!(int, pR)
	return WanIntijklR_pairindex!(int, (pair[1], pair[2], ipR))
end
function WanIntijklRHermitian_setindex!(int::WanIntijklRHermitian{T}, ::Val{true}, v::Number, ip1::Integer, ip2::Integer, iR::Integer) where {T}
	idx = findfirst(==((ip1, ip2, iR)), indices(int))
	isnothing(idx) || return int.vals[idx] = v
	push!(int.ips, ip1)
	push!(int.jps, ip2)
	push!(int.iRs, iR)
	push!(int.vals, v)
	ipR1 = int.pairpR[ip1]
	ipR2 = int.pairpR[ip2]
	conjR = int.R[iR] + int.pR[ipR2] - int.pR[ipR1]
	iconjR = get(int.R2idx, conjR, 0)
	if iszero(iconjR)
		push!(int.R, conjR)
		iconjR = int.nR
		int.R2idx[conjR] = iconjR
	end
	push!(int.iconjRs, iconjR)
	return v
end
function WanIntijklRHermitian_setindex!(int::WanIntijklRHermitian{T}, ::Val{false}, v::Number, ip1::Integer, ip2::Integer, iconjR::Integer) where {T}
	ip1 = int.conjpair[ip1]
	ip2 = int.conjpair[ip2]
	v = conj(v)
	ipR1 = int.pairpR[ip1]
	ipR2 = int.pairpR[ip2]
	R = int.R[iconjR] - (int.pR[ipR2] - int.pR[ipR1])
	iR = get(int.R2idx, R, 0)
	if iszero(iR)
		push!(int.R, R)
		iR = int.nR
		int.R2idx[R] = iR
	else
		idx = findfirst(==((ip1, ip2, iR)), indices(int))
		isnothing(idx) || return int.vals[idx] = v
	end
	push!(int.ips, ip1)
	push!(int.jps, ip2)
	push!(int.iRs, iR)
	push!(int.vals, v)
	push!(int.iconjRs, iconjR)
	return v
end
function WanIntijklR_setindex!(int::WanIntijklRHermitian{T}, v::Number, ip1::Integer, ip2::Integer, iR::Integer) where {T}
	p1t = int.pairtype[ip1]
	p2t = int.pairtype[ip2]
	return WanIntijklRHermitian_setindex!(int, Val((p1t, p2t) ∈ int.uplopair), v, ip1, ip2, iR)
end
