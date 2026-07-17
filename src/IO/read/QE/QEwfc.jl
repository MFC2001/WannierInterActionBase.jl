function Base.read(path::AbstractString, ::Type{QEwfc}; kwargs...)

	(_, ext) = splitext(path)

	if ext == ".dat"
		wfctype = QEwfc_dat
	elseif ext == ".hdf5"
		wfctype = QEwfc_hdf5
	else
		error("Unable to determine the format from the suffix of file.")
	end

	return read(path, wfctype; kwargs...)
end
function Base.read(io::IO, ::Type{QEwfc_dat}; convert = "native")

	io = FortranFile(io; convert)

	(ik, xk, ispin, gamma_only_int, scalef) = read(io, Int32, (Float64, 3), Int32, Int32, Float64)
	gamma_only = gamma_only_int ≠ 0

	(ngw, igwx, npol, nbnd) = read(io, (Int32, 4))

	rlattice = read(io, (Float64, 3, 3))

	xk = ReducedCoordinates{Float64}(rlattice \ xk)
	rlattice = ReciprocalLattice(rlattice / bohr_radius)

	G = map(ReducedCoordinates{Int}, eachcol(read(io, (Int32, 3, igwx))))

	n = npol * igwx
	components = Matrix{ComplexF64}(undef, n, nbnd)
	for i in 1:nbnd
		components[:, i] = read(io, (ComplexF64, n))
	end

	return QEwfc(ik, xk, ispin, gamma_only, scalef, ngw, igwx, npol, nbnd, rlattice, G, components)
end
function Base.read(path::AbstractString, ::Type{QEwfc_hdf5})

	wfc = h5open(path, "r")

	gamma_only = read_attribute(wfc, "gamma_only") == ".TRUE."
	igwx = read_attribute(wfc, "igwx")
	ik = read_attribute(wfc, "ik")
	ispin = read_attribute(wfc, "ispin")
	nbnd = read_attribute(wfc, "nbnd")
	ngw = read_attribute(wfc, "ngw")
	npol = read_attribute(wfc, "npol")
	scalef = read_attribute(wfc, "scale_factor")
	xk = read_attribute(wfc, "xk")

	G = map(ReducedCoordinates{Int}, eachcol(wfc["MillerIndices"][]))
	bg1 = read_attribute(wfc["MillerIndices"], "bg1")
	bg2 = read_attribute(wfc["MillerIndices"], "bg2")
	bg3 = read_attribute(wfc["MillerIndices"], "bg3")

	rlattice = [bg1 bg2 bg3]

	xk = ReducedCoordinates{Float64}(rlattice \ xk)
	rlattice = ReciprocalLattice(rlattice / bohr_radius)

	evc = wfc["evc"][]
	n = npol * igwx
	components = Matrix{ComplexF64}(undef, n, nbnd)
	for i in 1:nbnd, iG in 1:n
		iG2 = 2*iG
		components[iG, i] = complex(evc[iG2-1, i], evc[iG2, i])
	end

	close(wfc)

	return QEwfc(ik, xk, ispin, gamma_only, scalef, ngw, igwx, npol, nbnd, rlattice, G, components)
end
