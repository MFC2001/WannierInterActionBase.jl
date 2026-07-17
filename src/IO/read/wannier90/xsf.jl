function Base.read(io::IO, ::Type{wannier90_xsf}; period = [1, 1, 1], eps2::Real = 0)

	#comments
	for _ in 1:4
		readline(io)
	end

	#CRYSTAL
	readline(io)

	#PRIMVEC
	readline(io)
	lattice = Matrix{Float64}(undef, 3, 3)
	for i in 1:3
		lattice[:, i] = [parse(Float64, ss) for ss in split(readline(io))]
	end
	lattice = Lattice(lattice)

	#CONVVEC
	readline(io)
	convvec = Matrix{Float64}(undef, 3, 3)
	for i in 1:3
		convvec[:, i] = [parse(Float64, ss) for ss in split(readline(io))]
	end
	convvec = Lattice(convvec)

	#PRIMCOORD
	readline(io)
	(natom, _) = parse.(Int, split(readline(io)))
	atom_name = Vector{String}(undef, natom)
	atom_location = Vector{ReducedCoordinates{Float64}}(undef, natom)
	for i in 1:natom
		line = split(readline(io))
		atom_name[i] = line[1]
		atom_location_car = CartesianCoordinates{Float64}(parse.(Float64, line[2:4]))
		atom_location[i] = lattice \ atom_location_car
	end

	cell = Cell(lattice, atom_location, "Reduced"; name = atom_name, period)

	#BEGIN_BLOCK_DATAGRID_3D
	while !eof(io)
		line = readline(io)
		occursin(r"^\s*$", line) ? continue : break
	end
	#3D_field
	readline(io)
	#BEGIN_DATAGRID_3D_UNKNOWN
	readline(io)

	grid_size = NTuple{3, Int}(parse.(Int, split(readline(io))))
	origin = Vec3{Float64}(parse.(Float64, split(readline(io))))
	grid_vecs = Matrix{Float64}(undef, 3, 3)
	for i in 1:3
		grid_vecs[:, i] = [parse(Float64, ss) for ss in split(readline(io))]
	end
	grid_vecs = Mat3{Float64}(grid_vecs)

	value = Array{Float64, 3}(undef, grid_size...)
	nline = ceil(Int, prod(grid_size) // 6)
	for i in 1:(nline-1)
		value[((i-1)*6+1):(i*6)] = parse.(Float64, split(readline(io)))
	end
	value[((nline-1)*6+1):end] = parse.(Float64, split(readline(io)))

	#END_DATAGRID_3D
	readline(io)
	#END_BLOCK_DATAGRID_3D
	readline(io)

	# frame = SMatrix{3, 3}(frame)
	# dvec = SMatrix{3, 3}([frame[:, 1] / Nmesh[1] frame[:, 2] / Nmesh[2] frame[:, 3] / Nmesh[3]])

	# mesh = [repeat(0:Nmesh[1]-1, outer = Nmesh[2]) repeat(0:Nmesh[2]-1, inner = Nmesh[1])]
	# mesh = [repeat(mesh, outer = (Nmesh[3], 1)) repeat(0:Nmesh[3]-1, inner = Nmesh[1] * Nmesh[2])]

	# startpoint = round.(Int, inv(dvec) * startpoint)
	# mesh = transpose(mesh) .+ startpoint

	# dV = abs(dot(dvec[:, 3], cross(dvec[:, 1], dvec[:, 2])))
	# normconst = sum(abs2, value) * dV

	# I = abs2.(value) .> eps2
	# value = value[I]
	# mesh = mesh[:, I]

	return wannier90_xsf{Float64}(cell, origin, grid_vecs, grid_size, value)
end
