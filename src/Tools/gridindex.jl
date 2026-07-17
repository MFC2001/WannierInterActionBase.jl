"""
	gridindex(grid) -> Vector{ReducedCoordinates{Int}}

The examples of `grid` are [3,3,1], [-1:1, 0:0, -1:1], (3,3,3) or (-1:1,-1:1,-1:1).
Return a vector like [[0,0,0], [1,2,3],...].
"""
function gridindex(grid::AbstractVector{<:Integer})::Vector{ReducedCoordinates{Int}}

	length(grid) == 3 || throw(ArgumentError("should input 3 integers vector for grid"))

	start = -floor.(Int, (grid .- 1) .// 2)
	stop = ceil.(Int, (grid .- 1) .// 2)

	grid_index = [ReducedCoordinates(i, j, k) for i ∈ start[1]:stop[1], j ∈ start[2]:stop[2], k ∈ start[3]:stop[3]]

	return vec(grid_index)
end
function gridindex(grid::AbstractVector{<:UnitRange{<:Integer}})::Vector{ReducedCoordinates{Int}}

	length(grid) == 3 || throw(ArgumentError("should input 3 integers vector for grid"))

	all(x -> step(x) == 1, grid) || throw(ArgumentError("step should be 1 for each dimension in grid"))

	grid = UnitRange{Int}.(collect(grid))

	grid_index = [ReducedCoordinates(i, j, k) for i ∈ grid[1], j ∈ grid[2], k ∈ grid[3]]

	return vec(grid_index)
end
function gridindex(grid::NTuple{3, <:Integer})::Vector{ReducedCoordinates{Int}}

	start = .-floor.(Int, (grid .- 1) .// 2)
	stop = ceil.(Int, (grid .- 1) .// 2)

	grid_index = [ReducedCoordinates(i, j, k) for i ∈ start[1]:stop[1], j ∈ start[2]:stop[2], k ∈ start[3]:stop[3]]

	return vec(grid_index)
end
function gridindex(grid::NTuple{3, <:UnitRange{<:Integer}})::Vector{ReducedCoordinates{Int}}

	all(x -> step(x) == 1, grid) || throw(ArgumentError("step should be 1 for each dimension in grid"))

	grid = UnitRange{Int}.(grid)

	grid_index = [ReducedCoordinates(i, j, k) for i ∈ grid[1], j ∈ grid[2], k ∈ grid[3]]

	return vec(grid_index)
end
gridindex(sym::Symbol, args...; kwargs...) = gridindex(Val(sym), args...; kwargs...)
function gridindex(::Val{:WignerSeitz}, lattice, grid)
	points = MonkhorstPack(grid[1], grid[2], grid[3])
	points_inFBZ = fold2FBZ(lattice, points)
	points_inFBZ = map(points_inFBZ) do p
		ReducedCoordinates{Int}(p[1] * grid[1], p[2] * grid[2], p[3] * grid[3])
	end
	return points_inFBZ
end
"""
	gridindex(lattice, maxlattdist::Real, period::AbstractVector{Bool}) -> Vector{ReducedCoordinates{Int}}

Get all the lattice vectors whose norms are smaller than maxlattdist.
"""
function gridindex(lattice::Union{AbstractMatrix{<:Real}, AbstractLattice{<:Real}}, maxlattdist::Real, period::AbstractVector{Bool})::Vector{ReducedCoordinates{Int}}

	lattice = Lattice(parent(lattice))
	period = Vec3{Bool}(period)

	p = count(period)
	if p == 3
		ucpath = gridindex(lattice, maxlattdist, Val(3))
	elseif p == 2
		if !period[1]
			xyindex = SVector{2, Int}(2, 3)
		elseif !period[2]
			xyindex = SVector{2, Int}(3, 1)
		else
			xyindex = SVector{2, Int}(1, 2)
		end
		ucpath = gridindex(lattice, maxlattdist, xyindex, Val(2))
	elseif p == 1
		if period[1]
			zindex = 1
		elseif period[2]
			zindex = 2
		else
			zindex = 3
		end
		ucpath = gridindex(lattice, maxlattdist, zindex, Val(1))
	else
		ucpath = [ReducedCoordinates(0, 0, 0)]
	end
	return ucpath
end
function gridindex(lattice, maxlattdist, ::Val{3})

	𝐚 = lattice[begin:(begin+2)]
	𝐛 = lattice[(begin+3):(begin+5)]
	𝐜 = lattice[(begin+6):end]

	Ω = abs((𝐚 × 𝐛) ⋅ 𝐜)

	h₁ = Ω / norm(𝐛 × 𝐜)
	h₂ = Ω / norm(𝐜 × 𝐚)
	h₃ = Ω / norm(𝐚 × 𝐛)

	grid = gridindex(Int.(cld.(maxlattdist * 1.2, [h₁, h₂, h₃])) * 2 .+ 1)
	maxlattdist2 = maxlattdist^2
	ldot = parent(lattice)
	ldot = Mat3(transpose(ldot) * ldot)
	grid = filter(R -> R ⋅ (ldot * R) < maxlattdist2, grid)

	return grid
end
function gridindex(lattice, maxlattdist, xyindex, ::Val{2})

	𝐚 = lattice[begin:(begin+2)]
	𝐛 = lattice[(begin+3):(begin+5)]
	𝐜 = lattice[(begin+6):end]

	Ω = abs((𝐚 × 𝐛) ⋅ 𝐜)

	h₁ = Ω / norm(𝐛 × 𝐜)
	h₂ = Ω / norm(𝐜 × 𝐚)
	h₃ = Ω / norm(𝐚 × 𝐛)


	grid_3D = Int.(cld.(maxlattdist * 1.2, [h₁, h₂, h₃])) * 2 .+ 1
	grid = [1, 1, 1]
	grid[xyindex] .= grid_3D[xyindex]

	grid = gridindex(grid)
	maxlattdist2 = maxlattdist^2
	ldot = parent(lattice)
	ldot = Mat3(transpose(ldot) * ldot)
	grid = filter(R -> R ⋅ (ldot * R) < maxlattdist2, grid)

	return grid
end
function gridindex(lattice, maxlattdist, zindex, ::Val{1})

	𝐚 = lattice[begin:(begin+2)]
	𝐛 = lattice[(begin+3):(begin+5)]
	𝐜 = lattice[(begin+6):end]

	Ω = abs((𝐚 × 𝐛) ⋅ 𝐜)

	h₁ = Ω / norm(𝐛 × 𝐜)
	h₂ = Ω / norm(𝐜 × 𝐚)
	h₃ = Ω / norm(𝐚 × 𝐛)

	grid_3D = Int.(cld.(maxlattdist * 1.2, [h₁, h₂, h₃])) * 2 .+ 1
	grid = [1, 1, 1]
	grid[zindex] = grid_3D[zindex]

	grid = gridindex(grid)
	maxlattdist2 = maxlattdist^2
	ldot = parent(lattice)
	ldot = Mat3(transpose(ldot) * ldot)
	grid = filter(R -> R ⋅ (ldot * R) < maxlattdist2, grid)

	return grid
end
