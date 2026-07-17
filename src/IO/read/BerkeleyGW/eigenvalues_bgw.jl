function Base.read(io::IO, ::Type{eigenvalues_bgw})

	line = readline(io)
	neig = parse(Int, split(line, '=')[2])

	line = readline(io)
	vol = parse(Float64, split(line, '=')[2])

	line = readline(io)
	data = split(split(line, '=')[2])
	nspin = parse(Int, data[1])
	nspinor = parse(Int, data[2])

	readline(io)

	data = readdlm(io)

	eig = data[:, 1]
	dipole_abs2 = data[:, 2]
	dipole = complex.(data[:, 3], data[:, 4])
	return eigenvalues_bgw(neig, vol, nspin, nspinor, eig, dipole_abs2, dipole)
end
