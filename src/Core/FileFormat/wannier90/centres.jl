export wannier90_centres, numorb, spin_centres
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
	name = String[],
	index = Int[],
	atom_location = CartesianCoordinates{Float64}[],
	atom_name = String[],
	belonging = Int[],
)

	P = reduce(promote_type, eltype.(location))
	location = map(CartesianCoordinates{P} ∘ collect, location)

	P = reduce(promote_type, eltype.(atom_location))
	atom_location = map(CartesianCoordinates{P} ∘ collect, atom_location)

	num = length(location)

	name = deepcopy(name)
	n = length(name)
	if n ≠ 0 && n ≠ num
		error("Wrong orbital's name.")
	elseif eltype(name) <: AbstractString
		name = string.(name)
	end

	index = deepcopy(index)
	if isempty(index)
		index = collect(1:num)
	elseif length(index) ≠ num
		error("Wrong index of POSCAR.")
	end

	belonging = deepcopy(belonging)
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
		deepcopy(atom_name),
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
