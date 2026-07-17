"""
	findneighbors(cell::Cell, maxdist::Real; shell_tol::Real = 1e-6) -> (neighborlist, neighbordist)

Find neighbors in the form of shell.
The format of neighborlist: [ia1, ia2, Rx, Ry, Rz].
"""
function findneighbors(cell::Cell{P}, maxdist::Real; shell_tol::Real = 1e-6) where {P}

	h = maximum(1:3) do i
		a = cell.lattice[(i-1)*3+1:i*3]
		norm(a)
	end
	Rlist = gridindex(cell.lattice, maxdist + 4 * h, cell.period)
	Rlist_car = map(R -> cell.lattice * R, Rlist)

	if P <: ReducedCoordinates
		cell = convert(cell)
	end

	atomlocat_unit = cell.location

	natom = numatom(cell)
	nR = length(Rlist)

	atomlocat_super = vec([τ + R for τ in atomlocat_unit, R in Rlist_car])
	atomindex_super = vec([(ia, Rlist[iR]...) for ia in 1:natom, iR in 1:nR])

	# [na][nshell][na]
	neighborlist = Vector{Vector{Vector{NTuple{5, Int}}}}(undef, natom)
	neighbordist = Vector{Vector{Float64}}(undef, natom)

	for i in eachindex(atomlocat_unit)
		dist = [norm(l - atomlocat_unit[i]) for l in atomlocat_super]
		sort_idx = sortperm(dist)
		sort_idx = sort_idx[2:end]
		dist = dist[sort_idx]

		neighborlist_i = Vector{Vector{NTuple{5, Int}}}(undef, 1)
		neighborlist_i[1] = [(i, i, 0, 0, 0)]
		neighbordist_i = Vector{Float64}(undef, 1)
		neighbordist_i[1] = 0.0

		for (dist_i, dist_v) in zip(sort_idx, dist)
			if isapprox(dist_v, neighbordist_i[end]; atol = shell_tol)
				push!(neighborlist_i[end], (i, atomindex_super[dist_i]...))
			elseif dist_v > maxdist
				break
			elseif dist_v > neighbordist_i[end] + shell_tol
				push!(neighbordist_i, dist_v)
				push!(neighborlist_i, [(i, atomindex_super[dist_i]...)])
			else
				error("Something is wrong!")
			end
		end

		neighborlist[i] = neighborlist_i
		neighbordist[i] = neighbordist_i
	end

	return neighborlist, neighbordist
end
