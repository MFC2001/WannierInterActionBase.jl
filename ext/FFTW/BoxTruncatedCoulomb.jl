struct BoxTruncatedCoulomb <: TruncatedCoulomb
	lattice::Lattice{Float64}
	Glist::Vector{ReducedCoordinates{Int}}
	fft_grid::NTuple{3, Int}
	coulomb::Vector{Float64}
end
(TC::BoxTruncatedCoulomb)() = TC.coulomb
function TruncatedCoulomb(::Val{:box}, lattice::Union{Lattice, AbstractMatrix{<:Real}}, Glist;
	fft_grid, fft_nt::Integer = Threads.nthreads())

	check_lattice_period(lattice, [false, false, false])
	lattice = Lattice(parent(lattice))

	Ω = abs(det(parent(lattice)))

	half_lx = lattice[1, 1] / 2
	half_ly = lattice[2, 2] / 2
	half_lz = lattice[3, 3] / 2

	fft_grid = (Int(fft_grid[1]), Int(fft_grid[2]), Int(fft_grid[3]))
	coff = half_lx * half_ly * half_lz / (fft_grid[1] * fft_grid[2] * fft_grid[3]) * CoulombScale / Ω

	x = [half_lx * (i-0.5) / fft_grid[1] for i in 1:fft_grid[1]]
	y = [half_ly * (i-0.5) / fft_grid[2] for i in 1:fft_grid[2]]
	z = [half_lz * (i-0.5) / fft_grid[3] for i in 1:fft_grid[3]]
	invr = [1.0/sqrt(xx^2 + yy^2 + zz^2) for xx in x, yy in y, zz in z]

	# 后期可以将其调整为使用1Dfft+Julia层级并行。
	FFTW_nt = FFTW.get_num_threads()
	FFTW.set_num_threads(fft_nt)
	FFTW.r2r!(invr, FFTW.REDFT10)
	FFTW.set_num_threads(FFTW_nt)

	nG = length(Glist)
	coulomb = Vector{Float64}(undef, nG)
	for iG in Base.OneTo(nG)
		G = Glist[iG]
		G̃x, G̃y, G̃z = abs(G[1]) + 1, abs(G[2]) + 1, abs(G[3]) + 1
		coulomb[iG] = invr[G̃x, G̃y, G̃z] * coff
	end

	return BoxTruncatedCoulomb(lattice, Glist, fft_grid, coulomb)
end
