"""
	findneighbors(cell::Cell, maxdist::Real; shell_tol::Real = 1e-6) -> (neighborlist, neighbordist)

Find the neighbors of every atom in `cell` within the distance `maxdist`, grouped into shells of equal distance.

A neighbor of the atom `ia` is an atom `ja` in the unit cell translated by the lattice vector `R = (Rx, Ry, Rz)`, 
whose position is `τ[ja] + lattice * R`, where `τ` is the reduced coordinate of `ja`. 
Its distance to `ia` is therefore `norm(τ[ja] + lattice * R - τ[ia])`.

# Arguments
- `cell`: the unit cell. A cell storing [`ReducedCoordinates`](@ref) is converted to [`CartesianCoordinates`](@ref) internally.
- `maxdist`: the maximum distance to search, in the same length unit as `cell.lattice`.
- `shell_tol`: the absolute tolerance used to group distances into shells, i.e. distances differing by less than `shell_tol` belong to the same shell. Defaults to `1e-6`.

# Returns
- `neighborlist::Vector{Vector{Vector{NTuple{5, Int}}}}`: `neighborlist[ia][nsh]` collects the neighbors in the `nsh`-th shell of the atom `ia`. Each neighbor is a 5-tuple `(ia, ja, Rx, Ry, Rz)`.
- `neighbordist::Vector{Vector{Float64}}`: `neighbordist[ia][nsh]` is the distance of the `nsh`-th shell of the atom `ia`.

Both outer vectors have length `numatom(cell)`, and the shells are ordered by increasing distance, so `neighbordist[ia]` is ascending.

# Notes
- The first shell of every atom is the atom itself: `neighborlist[ia][1] == [(ia, ia, 0, 0, 0)]` with distance `0.0`. Drop `neighborlist[ia][1]` to obtain only the distinct neighbors.
- Every pair is reported from both sides: if `(ia, ja, R)` appears in `neighborlist[ia]`, then `(ja, ia, -R)` appears in `neighborlist[ja]`.
- For a non-periodic direction, the corresponding component of `R` is always zero.
- The candidate lattice vectors are enumerated with [`gridindex`](@ref) inside the sphere of radius `maxdist + 4h`, 
where `h` is the length of the longest basis vector. The cost therefore grows with `numatom(cell)` times the number of candidate vectors, so a large `maxdist` can be expensive.

# Examples
```julia
# A 2D kagome lattice whose nearest-neighbor distance is 1.
julia> cell = IdealCell(:kagome);
julia> neighborlist, neighbordist = findneighbors(cell, 1.0);

julia> neighbordist[1];     # [0.0, ≈1.0]: the self shell, then the nearest-neighbor shell
julia> neighborlist[1][2];  # the 4 nearest neighbors of atom 1, as (1, ja, Rx, Ry, Rz)
```
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
