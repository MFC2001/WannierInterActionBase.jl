function _blist_to_shell_weights(blist::AbstractVector; shell_tol::Real = 1e-6, orth_tol::Real = 1e-6)
	dist = map(bb->norm(bb), blist)
	dist_index = sortperm(dist)
	dist = dist[dist_index]

	dist_shell = Vector{Float64}(undef, 1)
	shell_bindex = Vector{Vector{Int}}(undef, 1)
	dist_shell[1] = dist[1]
	shell_bindex[1] = Int[]
	for (d_i, d_v) in zip(dist_index, dist)
		if isapprox(d_v, dist_shell[end]; atol = shell_tol)
			push!(shell_bindex[end], d_i)
		elseif d_v > dist_shell[end] + shell_tol
			push!(dist_shell, d_v)
			push!(shell_bindex, [d_i])
		else
			error("Something is wrong!")
		end
	end

	Amat = Matrix{Float64}(undef, 6, length(dist_shell))
	for (i, I) in enumerate(shell_bindex)
		Amat[:, i] = sum(blist[I]) do b
			[b[1]^2, b[2]^2, b[3]^2, b[1] * b[2], b[2] * b[3], b[3] * b[1]]
		end
	end
	F = svd(Amat)
	weights = F.V * Diagonal(1 ./ F.S) * F.U' * [1, 1, 1, 0, 0, 0]
	isapprox(Amat * weights, [1, 1, 1, 0, 0, 0]; atol = orth_tol) || error("Wrong b list!")

	weights_all = Vector{Float64}(undef, length(blist))
	for (i, w) in enumerate(weights)
		weights_all[shell_bindex[i]] .= w
	end

	return weights_all
end
