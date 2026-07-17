
export Cell, check_lattice_period, atomnames, numatom


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

Create a cell.

- `lattice`: a [`Lattice`](@ref) type or a matrix, note the basis vectors of the matrix are stored as columns.
- `location`: a vector contains all atom's location, its length equal to ``N``;
- `location_type`: \"Reduced\" or \"Cartesian\";
- `name`: a list of ``N`` values, where the same kind of name need to be a String.
- `index`: a list of ``N`` values, where the same kind of index need to be a Integer.
- `period`: [1(0), 1(0), 1(0)].

Use [`check_lattice_period`](@ref) to check whether the lattice match the period.
"""
function Cell(lattice, location::AbstractVector, location_type = "Cartesian";
	name::AbstractVector{<:Union{AbstractString, AbstractChar}} = String[],
	index::AbstractVector{<:Integer} = Int[],
	period = [true, true, true])

	check_lattice_period(lattice, period)

	period = Vec3{Bool}(period)

	num = length(location)

	n = length(name)
	if n ≠ 0 && n ≠ num
		throw(DimensionMismatch("The lengths of atomic locations and atomic names are different!"))
	end
	name = collect(string.(name))

	if isempty(index)
		index = collect(1:num)
	elseif length(index) ≠ num
		throw(DimensionMismatch("The lengths of atomic locations and atomic indexs are different!"))
	else
		index = collect(Int.(index))
	end

	lattice = Lattice(parent(lattice))

	P = reduce(promote_type, eltype.(location))
	if location_type[1] ∈ ['R', 'r']
		P = ReducedCoordinates{P}
	elseif location_type[1] ∈ ['C', 'c']
		P = CartesianCoordinates{P}
	else
		error("`location_type` only can be set as \"Reduced\" or \"Cartesian\".")
	end
	location = map(P ∘ collect, location)

	return Cell{P}(lattice, location, name, index, period)
end
"""
	check_lattice_period(lattice, period)

Check whether the lattice match the period.
For 2D caase, require the non-periodic direction should be perpendicular to the periodic plane;
For 1D and 0D case, require the basis vectors should be perpendicular to each other.
"""
function check_lattice_period(lattice, period)
	period = Vec3{Bool}(period)
	np = count(period)
	lattice = Mat3{Float64}(parent(lattice))
	if np == 2
		if !period[1]
			xyindex = SVector{2, Int}(2, 3)
			zindex = 1
		elseif !period[2]
			xyindex = SVector{2, Int}(3, 1)
			zindex = 2
		else
			xyindex = SVector{2, Int}(1, 2)
			zindex = 3
		end
		all(iszero, lattice[xyindex, zindex]) || error("Wrong lattice: the non-periodic direction should be perpendicular to the periodic plane.")
		all(iszero, lattice[zindex, xyindex]) || error("Wrong lattice: the non-periodic direction should be perpendicular to the periodic plane.")
	elseif np == 1 || np == 0
		is_diagonal(lattice) || error("Wrong lattice: the basis vectors should be perpendicular to each other.")
	end
	return nothing
end

atomnames(cell::Cell) = unique(cell.name)
numatom(cell::Cell) = length(cell.location)
Lattice(cell::Cell) = cell.lattice

function Base.sort(cell::Cell{P}) where {P}
	T = sortslices([cell.name cell.location collect(1:length(cell.name))]; dims = 1)

	name = similar(cell.name)
	name .= T[:, 1]

	location = similar(cell.location)
	location .= T[:, 2]

	index = cell.index[T[:, 3]]

	return Cell{P}(cell.lattice, location, name, index, cell.period)
end

function Base.convert(cell::Cell{ReducedCoordinates{P}}) where {P}

	location = map(x -> cell.lattice * x, cell.location)

	return Cell{CartesianCoordinates{Float64}}(cell.lattice, location, copy(cell.name), copy(cell.index), cell.period)
end
function Base.convert(cell::Cell{CartesianCoordinates{P}}) where {P}

	location = map(x -> cell.lattice \ x, cell.location)

	return Cell{ReducedCoordinates{Float64}}(cell.lattice, location, copy(cell.name), copy(cell.index), cell.period)
end
function Base.convert(::Type{Cell{P}}, cell::Cell) where {P}
	return Cell{P}(cell.lattice, P.(location), copy(cell.name), copy(cell.index), cell.period)
end
function Base.zero(::Type{Cell{P}}) where {P}
	return Cell{P}(one(Lattice{Float64}), P[], String[], Int[], Bool[1, 1, 1])
end
function Base.zero(::Type{Cell})
	return Cell{ReducedCoordinates{Int}}(one(Lattice{Float64}), ReducedCoordinates{Int}[], String[], Int[], Bool[1, 1, 1])
end
