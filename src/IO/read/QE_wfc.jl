
function Base.read(io::IO, ::Type{QE_wfc}; convert = "native")

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

	return QE_wfc(ik, xk, ispin, gamma_only, scalef, ngw, igwx, npol, nbnd, rlattice, G, components)
end
