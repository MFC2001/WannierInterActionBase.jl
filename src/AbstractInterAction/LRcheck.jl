export LRcheck
"""
	LRcheck(V, U::HR, lattice::Lattice, orblocat::AbstractVector{<:AbstractVector{<:Real}}, longrcut::Real)

Obtain the maximum deviation of U-type interactions beyond longrcut and V(r).
Note: V(r) = V(norm(lattice * r_frac)), r_frac = R + τⱼ - τᵢ.
"""
function LRcheck(V, U::HR, lattice::Lattice, orblocat::AbstractVector{<:AbstractVector{<:Real}}, longrcut::Real)

	lattice = parent(lattice)
	ldot = Mat3(transpose(lattice) * lattice)

	longrcut2 = longrcut * longrcut
	(maxeps, maxeps_idx) = findmax(Base.OneTo(length(U.value))) do i
		path = ReducedCoordinates(U.path[1, i], U.path[2, i], U.path[3, i])
		i_idx = U.path[4, i]
		j_idx = U.path[5, i]
		r_frac = path + orblocat[j_idx] - orblocat[i_idx]
		r2 = r_frac ⋅ (ldot * r_frac)
		r2 ≤ longrcut2 ? 0.0 : abs(V(sqrt(r2)) - U.value[i])
	end
	path = ReducedCoordinates(U.path[1, maxeps_idx], U.path[2, maxeps_idx], U.path[3, maxeps_idx])
	i_idx = U.path[4, maxeps_idx]
	j_idx = U.path[5, maxeps_idx]

	return (path, i_idx, j_idx, maxeps)
end
