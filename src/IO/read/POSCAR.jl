
function Base.read(io::IO, ::Type{POSCAR}; period = Bool[1, 1, 1])::Cell
	#first line is comment.
	readline(io)

	scale = parse(Float64, readline(io))
	lattice = Matrix{Float64}(undef, 3, 3)
	for i in 1:3
		lattice[i, :] = [parse(Float64, ss) for ss in split(readline(io))]
	end

	elem_types = split(readline(io))
	elem_nums = [parse(Int, nn) for nn in split(readline(io))]

	str = readline(io)
	if str[1] ∈ ['S', 's'] #Slective Dynamics
		str = readline(io)
	end

	if str[1] ∈ ['D', 'd'] #Direct
		location_type = "Reduced"
	elseif str[1] ∈ ['C', 'c'] #Cartesian
		location_type = "Cartesian"
	else
		error("Wrong POSCARfile.")
	end

	natoms = sum(elem_nums)
	location = Vector{Vec3{Float64}}(undef, natoms)
	for i in eachindex(location)
		location[i] = Vec3(parse.(Float64, split(readline(io))[1:3]))
	end

	lattice = Lattice(scale * transpose(lattice))

	name = String[]
	for (elem_type, elem_num) in zip(elem_types, elem_nums)
		append!(name, fill(elem_type, elem_num))
	end

	return Cell(lattice, location, location_type; name, period)
end
