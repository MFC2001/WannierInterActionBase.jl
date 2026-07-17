"""
	fold2FBZ(lattice, points::AbstractVector{<:AbstractVector}) -> AbstractVector{<:AbstractVector}
	fold2FBZ(lattice, point::AbstractVector{<:Real}) -> AbstractVector{<:Real}

input reduced coordinates.
Return the input's equivalent points in FBZ.
"""
function fold2FBZ(lattice::Union{AbstractLattice, AbstractMatrix}, points::AbstractVector{V}; tol::Real = 1e-6) where {T <: Real, V <: AbstractVector{T}}

	FBZ = WignerSeitz(lattice)

	foldedkpoints = deepcopy(points)
	for index in eachindex(foldedkpoints)

		point = foldedkpoints[index]
		if _is_approx_integer(point; atol = 1e-6)
			foldedkpoints[index] = point - round.(Int, point)
			continue
		end

		point_coord = lattice * point

		while true
			(inFBZ, nΓ) = is_inside_detail(point_coord, FBZ)
			if inFBZ
				break
			else
				point_coord = point_coord - nΓ
			end
		end

		point = lattice \ point_coord
		if T <: Rational
			point = rationalize.(point; tol)
		end
		foldedkpoints[index] = point
	end

	return foldedkpoints
end

function fold2FBZ(lattice::Union{AbstractLattice, AbstractMatrix}, point::V; tol::Real = 1e-6) where {T <: Real, V <: AbstractVector{T}}

	FBZ = WignerSeitz(lattice)

	if _is_approx_integer(point; atol = 1e-6)
		return point - round.(point)
	end

	point_coord = lattice * point

	while true
		(inFBZ, nΓ) = is_inside_detail(point_coord, FBZ)
		if inFBZ
			break
		else
			point_coord = point_coord - nΓ
		end
	end

	point = lattice \ point_coord
	if T <: Rational
		point = rationalize.(point; tol)
	end

	return V(point)
end

fold2FBZ(lattice, kpoints::AbstractBrillouinZone) =
	fold2FBZ(lattice, kpoints.kdirect)
fold2FBZ(lattice, kpoints::MonkhorstPack) =
	fold2FBZ(lattice, RedKgrid(kpoints).kdirect)
