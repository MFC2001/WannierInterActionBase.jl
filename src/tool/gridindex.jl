export gridindex
"""
	gridindex(grid::AbstractVector{<:Integer}; n::Integer = 1)

	Return a grid, its datatype is like [[0,0,0], [1,2,3],...].
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
