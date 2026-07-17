export wannier90_centres, numorb, spin_centres, translate, translate!
"""
wannier90_centres.dat

	read(path/io, ::Type{wannier90_centres}) -> wannier90_centres

Read `wannier90_centres.dat`.

	write(path/io, orbital::wannier90_centres, ::Type{wannier90_centres})

Write `wannier90_centres.dat`.

	wannier90_centres <: FileFormat

The data in wannier90_centres.xyz. Fields:
- `location::Vector{CartesianCoordinates{Float64}}`: location of orbital centres;
- `name::Vector{String}`: name of orbitals;
- `index::Vector{Int}`: index of orbitals;
- `atom_location::Vector{CartesianCoordinates{Float64}}`: location of atoms;
- `atom_name::Vector{String}`: name of atoms;
- `belonging::Vector{Int}`: the index of atom that orbital belongs to.
"""
struct wannier90_centres <: FileFormat
	location::Vector{CartesianCoordinates{Float64}}
	name::Vector{String}
	index::Vector{Int}
	atom_location::Vector{CartesianCoordinates{Float64}}
	atom_name::Vector{String}
	belonging::Vector{Int} # The orbitals is belongs to which atom.
end
"""
	wannier90_centres(location::AbstractVector; name = String[], index = Int[], 
		atom_location = CartesianCoordinates{Float64}[], atom_name = String[], belonging = Int[])

Construct an `wannier90_centres` object.

- `location`: a vector contains all the locations of orbitals, its length equals ``N``;
- `name` and `index`: its element is the name or index of each orbital;
- `atom_location` and `atom_name`: its element is the location or name of each orbital;
- `belonging`: the i-th orbital belongs to the `belonging[i]`-th atom.
"""
function wannier90_centres(
	location::AbstractVector;
	name::AbstractVector{<:Union{AbstractString, AbstractChar}} = String[],
	index::AbstractVector{<:Integer} = Int[],
	atom_location::AbstractVector = CartesianCoordinates{Float64}[],
	atom_name::AbstractVector{<:Union{AbstractString, AbstractChar}} = String[],
	belonging::AbstractVector{<:Integer} = Int[],
)

	P = reduce(promote_type, eltype.(location))
	location = map(CartesianCoordinates{P} ∘ collect, location)

	P = reduce(promote_type, eltype.(atom_location))
	atom_location = map(CartesianCoordinates{P} ∘ collect, atom_location)

	num = length(location)

	n = length(name)
	if n ≠ 0 && n ≠ num
		throw(DimensionMismatch("The lengths of orbitals' locations and orbitals' names are different!"))
	end
	name = collect(string.(name))

	if isempty(index)
		index = collect(1:num)
	elseif length(index) ≠ num
		throw(DimensionMismatch("The lengths of orbitals' locations and orbitals' indices are different!"))
	else
		index = collect(Int.(index))
	end

	n = length(atom_name)
	if n ≠ 0 && n ≠ length(atom_location)
		throw(DimensionMismatch("The lengths of atomic locations and atomic names are different!"))
	end
	atom_name = collect(string.(atom_name))

	belonging = copy(belonging)
	if isempty(belonging) && !isempty(atom_location)
		belonging = Vector{Int}(undef, num)
		for i in 1:num
			(_, I) = findmin(x -> sum(abs2, x - location[i]), atom_location)
			belonging[i] = I
		end
	end

	return wannier90_centres(
		location,
		name,
		index,
		atom_location,
		atom_name,
		belonging,
	)
end

function Base.show(io::IO, orbital::wannier90_centres)
	print(io, "wannier90_centres with $(numorb(orbital)) orbitals and $(length(orbital.atom_location)) atoms.")
end
numorb(orbital::wannier90_centres) = length(orbital.location)

function spin_centres(orbital::wannier90_centres)::wannier90_centres

	belonging = orbital.belonging
	if !isempty(belonging)
		spinbelonging = [belonging; belonging]
	end

	return wannier90_centres(
		[orbital.location; orbital.location],
		[orbital.name; orbital.name],
		[orbital.index; orbital.index],
		orbital.atom_location,
		orbital.atom_name,
		spinbelonging,
	)
end

"""
	translate(orbital::wannier90_centres, lattice)

Translate orbital to [0, 0, 0] cell.

See also `translate!`.
"""
function translate(orbital::wannier90_centres, lattice)
	lattice = Lattice(parent(lattice))
	orblocation_frac = map(x -> lattice \ x, orbital.location)
	frac_incell = map(x -> mod.(x, 1), orblocation_frac)
	centres = map(eachindex(frac_incell)) do i
		centre = round.(Int, orblocation_frac[i] - frac_incell[i])
		return i => centre
	end
	return translate(orbital, centres...), centres
end
"""
	translate!(orbital::wannier90_centres, lattice)

Translate orbital to [0, 0, 0] cell.

See also `translate!`.
"""
function translate!(orbital::wannier90_centres, lattice)
	lattice = Lattice(parent(lattice))
	orblocation_frac = map(x -> lattice \ x, orbital.location)
	frac_incell = map(x -> mod.(x, 1), orblocation_frac)
	centres = map(eachindex(frac_incell)) do i
		centre = round.(Int, orblocation_frac[i] - frac_incell[i])
		return i => centre
	end
	return translate!(orbital, centres...), centres
end
"""
	translate(orbital::wannier90_centres, lattice, centres::Pair{<:Integer, <:AbstractVector{<:Integer}}...) -> wannier90_centres

You only need to input the wannier centres R which isn't in the [0, 0, 0] cell.

See also `translate!`.
"""
function translate(orbital::wannier90_centres, lattice,
	centres::Pair{<:Integer, <:AbstractVector{<:Integer}}...)::wannier90_centres
	lattice = Lattice(parent(lattice))
	orbital = deepcopy(orbital)
	for (iw, R) in centres
		orbital.location[iw] = orbital.location[iw] - lattice * R
	end
	return orbital
end
"""
	translate!(orbital::wannier90_centres, lattice, centres::Pair{<:Integer, <:AbstractVector{<:Integer}}...) -> wannier90_centres

You only need to input the wannier centres R which isn't in the [0, 0, 0] cell.

See also `translate`.
"""
function translate!(orbital::wannier90_centres, lattice,
	centres::Pair{<:Integer, <:AbstractVector{<:Integer}}...)::wannier90_centres
	lattice = Lattice(parent(lattice))
	for (iw, R) in centres
		orbital.location[iw] = orbital.location[iw] - lattice * R
	end
	return orbital
end
