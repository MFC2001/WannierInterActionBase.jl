function WanIntijklRHermitian_pre_pair!(ips, jps, pair, pR)
	# obtain pR and conjpR, and original_pR_idx
	original_pR = pR
	conj_pR = [-RR for RR in original_pR]
	pR = sort(unique([original_pR; conj_pR]); by = norm2)
	pR2idx = Dict(RR=>i for (i, RR) in enumerate(pR))
	original_pR_idx = [pR2idx[RR] for RR in original_pR]
	conjpR = [pR2idx[-RR] for RR in pR]
	# obtain pair and conjpair
	original_pair = [(p[1], p[2], original_pR_idx[p[3]]) for p in pair]
	conj_pair = [(p[2], p[1], conjpR[p[3]]) for p in original_pair]
	by(p) = (p[3], p[1], p[2])
	pair = sort(unique([original_pair; conj_pair]); by)
	pair2idx = Dict(p=>i for (i, p) in enumerate(pair))
	original_pair_idx = [pair2idx[p] for p in original_pair]
	conjpair = [pair2idx[(p[2], p[1], conjpR[p[3]])] for p in pair]
	# update
	@inbounds @simd for i in eachindex(ips)
		ips[i] = original_pair_idx[ips[i]]
	end
	@inbounds @simd for i in eachindex(jps)
		jps[i] = original_pair_idx[jps[i]]
	end
	return ips, jps, pair, pair2idx, conjpair, pR, pR2idx, conjpR
end
function WanIntijklRHermitian_pre_iconjR!(ips, jps, iRs, pair, pR, R, pairtype, uplopair)
	# R′
	pR_minus_pR = [R₁ - R₂ for R₁ in pR, R₂ in pR]
	R′ = sort(unique(vec(pR_minus_pR)); by = norm2)
	R′2idx = Dict(RR => i for (i, RR) in enumerate(R′))
	pR_minus_pR_to_R′ = map(RR -> R′2idx[RR], pR_minus_pR)
	pair_minus_pair_to_R′ = [pR_minus_pR_to_R′[p1[3], p2[3]] for p1 in pair, p2 in pair]
	# obtain R
	original_R = R
	conj_R = unique(vec([RR+RR′ for RR in original_R, RR′ in R′]))
	R = sort(unique([original_R; conj_R]); by = norm2)
	R2idx = Dict(RR => i for (i, RR) in enumerate(R))
	original_R_idx = [R2idx[RR] for RR in original_R]
	R_add_R′_to_R = [get(R2idx, RR + RR′, 0) for RR in R, RR′ in R′]
	# update
	iconjRs = Vector{Int}(undef, length(iRs))
	@inbounds @simd for i in eachindex(iRs)
		ip1 = ips[i]
		ip2 = jps[i]
		(pairtype[ip1], pairtype[ip2]) ∈ uplopair || error("Found pair interaction out of uplopair!")
		iR = original_R_idx[iRs[i]]
		iR′ = pair_minus_pair_to_R′[ip2, ip1]
		iRs[i] = iR
		iconjRs[i] = R_add_R′_to_R[iR, iR′]
	end
	return iRs, iconjRs, R, R2idx
end
function WanIntijklRHermitian_pre_value(ips, jps, iRs, vals::AbstractVector{T}, pR, R, conjpair, pairtype, pairpR, uplopair) where {T}
	# obtain R and all value.
	R2idx = Dict(RR => i for (i, RR) in enumerate(R))
	all_value = Dict{NTuple{3, Int}, Tuple{T, Int}}()
	@inbounds for (ip1, ip2, iR, v) in zip(ips, jps, iRs, vals)
		if (pairtype[ip1], pairtype[ip2]) ∈ uplopair
			ppR = (ip1, ip2, iR)
		else
			ip1 = conjpair[ip1]
			ip2 = conjpair[ip2]
			ipR1 = pairpR[ip1]
			ipR2 = pairpR[ip2]
			R_v = R[iR] - (pR[ipR2] - pR[ipR1])
			iR = get(R2idx, R_v, 0)
			if iszero(iR)
				push!(R, R_v)
				iR = length(R)
				R2idx[R_v] = iR
			end
			ppR = (ip1, ip2, iR)
			v = conj(v)
		end
		(v_sum, nv) = get(all_value, ppR, (zero(T), 0))
		all_value[ppR] = (v_sum + v, nv + 1)
	end
	# construct new data.
	n = length(all_value)
	new_ips = Vector{Int}(undef, n)
	new_jps = Vector{Int}(undef, n)
	new_iRs = Vector{Int}(undef, n)
	new_vals = Vector{ComplexF64}(undef, n)
	for (idx, (key, (v_sum, nv))) in enumerate(all_value)
		new_ips[idx] = key[1]
		new_jps[idx] = key[2]
		new_iRs[idx] = key[3]
		new_vals[idx] = v_sum / nv
	end
	# update Ridx.
	R_exist = sort(unique(new_iRs))
	new_R = R[R_exist]
	@inbounds if length(new_R) < length(R)
		R_map = Dict(iR=>i for (i, iR) in enumerate(R_exist))
		R_map = [get(R_map, iR, 0) for iR in eachindex(R)]
		@simd for i in eachindex(new_iRs)
			new_iRs[i] = R_map[new_iRs[i]]
		end
	end
	return new_ips, new_jps, new_iRs, new_vals, new_R
end
