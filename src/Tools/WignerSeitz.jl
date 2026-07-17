function WignerSeitz(lattice::AbstractLattice)
	neighbour = map(_neighbour_offsets_3d) do Δ
		CartesianCoordinates(lattice * Δ)
	end
	offsets = map(car -> 0.5 * (car[1] * car[1] + car[2] * car[2] + car[3] * car[3]), neighbour)
	WS = HalfspaceConvexPolyhedron(neighbour, offsets)
	return WS
end
# function WignerSeitz(lattice::AbstractLattice, period::AbstractVector{<:Integer})
# 	length(period) == 3 || error("Wrong length of period!")
# 	if count(period) == 3
# 		WignerSeitz(lattice)
# 	elseif count(period) == 2
# 		_WignerSeitz_2D(lattice, period)
# 	elseif count(period) == 1
# 		_WignerSeitz_1D(lattice, period)
# 	end
# 	return WS
# end
# function _WignerSeitz_2D(lattice::AbstractLattice, period::AbstractVector{<:Integer})
# 	if !period[1]
# 		xyindex = SVector{2, Int}(2, 3)
# 	elseif !period[2]
# 		xyindex = SVector{2, Int}(3, 1)
# 	else
# 		xyindex = SVector{2, Int}(1, 2)
# 	end
# 	lattice = lattice[xyindex, xyindex]
# 	return _WignerSeitz_2D(lattice)
# end
# function _WignerSeitz_1D(lattice::AbstractLattice, period::AbstractVector{<:Integer})
# 	if period[1]
# 		lattice = lattice[1, 1]
# 	elseif !period[2]
# 		lattice = lattice[2, 2]
# 	else
# 		lattice = lattice[3, 3]
# 	end
# 	return _WignerSeitz_1D(lattice)
# end

function WignerSeitz(lattice::AbstractMatrix{T}) where {T <: Real}
	if size(lattice) == (3, 3)
		return WignerSeitz(Lattice(lattice))
	elseif size(lattice) == (2, 2)
		return _WignerSeitz_2D(lattice)
	elseif size(lattice) == (1, 1)
		return _WignerSeitz_1D(lattice[1])
	else
		error("Only support 2D and 3D.")
	end
	return WS
end
function _WignerSeitz_2D(lattice::AbstractMatrix{T}) where {T <: Real}
	lattice = SMatrix{2, 2, T, 4}(lattice)
	neighbour = map(_neighbour_offsets_2d) do Δ
		SVector{2, T}(lattice * Δ)
	end
	offsets = map(car -> 0.5 * (car[1] * car[1] + car[2] * car[2]), neighbour)
	WS = HalfspaceConvexPolyhedron(neighbour, offsets)
	return WS
end
function _WignerSeitz_1D(lattice::T) where {T <: Real}
	error("Maybe you can do it by yourself!")
end



# Lattice方法歧义太多，仅对应3D。
