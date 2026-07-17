function Base.read(io::IO, ::Type{wannier90_nnkp})

	# real_lattice
	while true
		if readline(io) == "begin real_lattice"
			break
		end
	end
	lattice = Matrix{Float64}(undef, 3, 3)
	for i in 1:3
		lattice[:, i] = [parse(Float64, ss) for ss in split(readline(io))]
	end
	lattice = Lattice(lattice)

	# recip_lattice
	while true
		if readline(io) == "begin recip_lattice"
			break
		end
	end
	rlattice = Matrix{Float64}(undef, 3, 3)
	for i in 1:3
		rlattice[:, i] = [parse(Float64, ss) for ss in split(readline(io))]
	end
	rlattice = ReciprocalLattice(rlattice)

	# kpoints
	while true
		if readline(io) == "begin kpoints"
			break
		end
	end
	nk = parse(Int, readline(io))
	kpoints = Vector{ReducedCoordinates{Float64}}(undef, nk)
	for ik in 1:nk
		kpoints[ik] = ReducedCoordinates(parse.(Float64, split(readline(io))))
	end

	# nnkpts
	while true
		if readline(io) == "begin nnkpts"
			break
		end
	end
	nntot = parse(Int, readline(io))
	nnkpts = Matrix{Int}(undef, 5, nk*nntot)
	for ikb in 1:(nk*nntot)
		nnkpts[:, ikb] = parse.(Int, split(readline(io)))
	end

	# exclude_bands
	while true
		if readline(io) == "begin exclude_bands"
			break
		end
	end
	neb = parse(Int, readline(io))
	exclude_bands = Vector{Int}(undef, neb)
	for ib in 1:neb
		exclude_bands[ib] = parse(Int, readline(io))
	end

	return wannier90_nnkp(lattice, rlattice, kpoints, nnkpts, exclude_bands)
end
