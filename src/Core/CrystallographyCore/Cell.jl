
export Cell, atomnames, numatom


abstract type AbstractCell end

"""
	Cell{P <: Union{CartesianCoordinates{<:Real}, ReducedCoordinates{<:Real}}}

Contains all the information of a crystalline unitcell. Fields:
- `lattice::Lattice{Float64}`: lattice basis matrix;
- `location::Vector{P}`: all atoms' location;
- `name::Vector{String}`: all atoms' location;
- `index::Vector{Int}`: all atoms' index;
- `period::Vec3{Bool}`: the periodicty of unitcell.
"""
@struct_hash_equal_isequal struct Cell{P <: Union{CartesianCoordinates{<:Real}, ReducedCoordinates{<:Real}}} <: AbstractCell
	lattice::Lattice{Float64}
	location::Vector{P}
	name::Vector{String}
	index::Vector{Int}
	period::Vec3{Bool}
end

"""
	Cell(lattice, location, location_type = "Cartesian"; name = String[], index = Int[], period = Bool[1, 1, 1])

Create a new cell.

- `lattice`: a [`Lattice`](@ref) type or a matrix, note the basis vectors of the matrix are stored as columns.
- `location`: a vector contains all atom's location, its length equal to ``N``;
- `location_type`: \"Reduced\" or \"Cartesian\";
- `name`: a list of ``N`` values, where the same kind of name need to be a String.
- `index`: a list of ``N`` values, where the same kind of index need to be a Integer.
- `period`: [1(0), 1(0), 1(0)].

Make sure the basis at unperiodic direction is perpendicular to the other basis.
"""
function Cell(lattice, location::AbstractVector, location_type = "Cartesian"; name = String[], index = Int[], period = [1, 1, 1])

	num = length(location)

	n = length(name)
	if n ≠ 0 && n ≠ num
		throw(DimensionMismatch("The lengths of atomic locations and atomic names are different!"))
	end
	name = string.(name)

	index = deepcopy(index)
	if isempty(index)
		index = collect(1:num)
	elseif length(index) ≠ num
		throw(DimensionMismatch("The lengths of atomic locations and atomic indexs are different!"))
	end

	lattice = deepcopy(lattice)
	if !(lattice isa Lattice)
		lattice = Lattice(lattice)
	end


	P = reduce(promote_type, eltype.(location))
	if location_type[1] ∈ ['R', 'r']
		P = ReducedCoordinates{P}
	elseif location_type[1] ∈ ['C', 'c']
		P = CartesianCoordinates{P}
	else
		error("`location_type` only can be set as \"Reduced\" or \"Cartesian\".")
	end
	location = map(P ∘ collect, location)

	period = Vec3{Bool}(collect(period))

	return Cell{P}(lattice, location, name, index, period)
end

atomnames(cell::Cell) = unique(cell.name)
numatom(cell::Cell) = length(cell.location)
Lattice(cell::Cell) = deepcopy(cell.lattice)

function Base.sort(cell::Cell{P}) where {P}
	T = sortslices([cell.name cell.location collect(1:length(cell.name))]; dims = 1)

	name = similar(cell.name)
	name .= T[:, 1]

	location = similar(cell.location)
	location .= T[:, 2]

	index = cell.index[T[:, 3]]

	return Cell{P}(deepcopy(cell.lattice), location, name, index, deepcopy(cell.period))
end

function Base.convert(cell::Cell{ReducedCoordinates{P}}) where {P <: Real}

	location = map(x -> cell.lattice * x, cell.location)

	return Cell{CartesianCoordinates{P}}(deepcopy(cell.lattice), location, deepcopy(cell.name), deepcopy(cell.index), deepcopy(cell.period))
end
function Base.convert(cell::Cell{CartesianCoordinates{P}}) where {P <: Real}

	location = map(x -> cell.lattice \ x, cell.location)

	return Cell{ReducedCoordinates{P}}(deepcopy(cell.lattice), location, deepcopy(cell.name), deepcopy(cell.index), deepcopy(cell.period))
end
