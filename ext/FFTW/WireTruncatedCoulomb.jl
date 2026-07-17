struct WireTruncatedCoulomb{dctT <: FFTW.AbstractFFTs.Plan{Float64}} <: TruncatedCoulomb
	lattice::Lattice{Float64}
	period::Vec3{Bool}
	Glist::Vector{ReducedCoordinates{Int}}
	fft_grid::NTuple{2, Int}
	zindex::Int
	xyindex::SVector{2, Int}
	rlatticez::Float64
	Gz::Vector{Int}
	Gxy::Vector{SVector{2, Int}}
	Gzcar::Vector{Float64}
	dct_K0::dctT
	r::Matrix{Float64}
	coff::Float64
end
function (TC::WireTruncatedCoulomb)(q::AbstractVector{<:Real}, Gidx::AbstractVector{<:Integer})
	coulomb = Vector{Float64}(undef, length(Gidx))
	return TC(coulomb, q, Gidx)
end
@inbounds function (TC::WireTruncatedCoulomb)(coulomb::AbstractVector{<:Real}, q::AbstractVector{<:Real}, Gidx::AbstractVector{<:Integer})
	qcar = TC.rlatticez * q[TC.zindex]
	Gz_q = view(TC.Gz, Gidx)
	Gxy_q = view(TC.Gxy, Gidx)
	Gzcar_q = view(TC.Gzcar, Gidx)

	# unique Gz, and its index in Glist_q
	Gz_unique, Gz_indices = unique_groups(Gz_q)

	Threads.@threads for iGz in eachindex(Gz_unique)
		realgrid = Matrix{Float64}(undef, TC.fft_grid[1], TC.fft_grid[2])
		iGs_iGz = Gz_indices[iGz]
		qGz = qcar + Gzcar_q[iGs_iGz[1]]
		realgrid .= besselk.(0, qGz .* TC.r)
		TC.dct_K0 * realgrid
		for iG in iGs_iGz
			Gxy = Gxy_q[iG]
			G̃x, G̃y = abs(Gxy[1]) + 1, abs(Gxy[2]) + 1
			coulomb[iG] = realgrid[G̃x, G̃y] * TC.coff
		end
	end

	# iGz_chunks = _split_tasks(length(Gz_unique), TC.nchunk)
	# Threads.@threads for ichunk in Base.OneTo(TC.nchunk)
	# 	realgrid = TC.realgrids[ichunk]
	# 	for iGz in iGz_chunks[ichunk]
	# 		iGs_iGz = Gz_indices[iGz]
	# 		qGz = qcar + Gzcar_q[iGs_iGz[1]]
	# 		realgrid .= besselk.(0, qGz .* TC.r)
	# 		TC.dct_K0 * realgrid
	# 		qGz2 = qGz * qGz
	# 		for iG in iGs_iGz
	# 			Gxy = Gxy_q[iG]
	# 			G̃x, G̃y = abs(Gxy[1]) + 1, abs(Gxy[2]) + 1
	# 			coulomb[iG] = realgrid[G̃x, G̃y] * TC.coff
	# 		end
	# 	end
	# end

	return coulomb
end
function TruncatedCoulomb(::Val{:wire}, lattice::Union{Lattice, AbstractMatrix{<:Real}}, Glist;
	period = [false, false, true], fft_grid)

	check_lattice_period(lattice, period)
	period = Vec3{Bool}(period)
	@assert count(period) == 1 "Wrong period, wire truncation is used in 1D system."
	if period[1]
		zindex = 1
		xyindex = SVector{2, Int}(2, 3)
	elseif period[2]
		zindex = 2
		xyindex = SVector{2, Int}(3, 1)
	else
		zindex = 3
		xyindex = SVector{2, Int}(1, 2)
	end

	lattice = Lattice(parent(lattice))
	rlattice = reciprocal(lattice)

	Ω = abs(det(parent(lattice)))

	rlatticez = rlattice[zindex, zindex]

	Gz = map(Glist) do G
		G[zindex]
	end
	Gxy = map(Glist) do G
		G[xyindex]
	end
	Gzcar = Gz .* rlatticez

	half_lx = lattice[xyindex[1], xyindex[1]] / 2
	half_ly = lattice[xyindex[2], xyindex[2]] / 2

	fft_grid = (Int(fft_grid[xyindex[1]]), Int(fft_grid[xyindex[2]]))
	x = [half_lx * (i-0.5) / fft_grid[1] for i in 1:fft_grid[1]]
	y = [half_ly * (i-0.5) / fft_grid[2] for i in 1:fft_grid[2]]
	r = [sqrt(xx^2 + yy^2) for xx in x, yy in y]

	coff = half_lx * half_ly / (fft_grid[1] * fft_grid[2]) * 2 * CoulombScale / Ω

	FFTW_nt = FFTW.get_num_threads()
	FFTW.set_num_threads(1)
	# plan会保持建立时的线程数运行。
	dct_K0 = FFTW.plan_r2r!(Matrix{Float64}(undef, fft_grid[1], fft_grid[2]), FFTW.REDFT10)
	FFTW.set_num_threads(FFTW_nt)

	return WireTruncatedCoulomb(lattice, period, Glist, fft_grid,
		zindex, xyindex, rlatticez, Gz, Gxy, Gzcar, dct_K0, r, coff)
end
