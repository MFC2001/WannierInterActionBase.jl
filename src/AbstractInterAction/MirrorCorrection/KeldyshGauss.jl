# 适用于2D系统。
function MirrorCorrection(::Val{:keldysh_gauss_optim}, U::AbstractReciprocalHoppings{T}, kgrid::RedKgrid,
	lattice::Lattice, orblocat::Vector{<:ReducedCoordinates}; norb::Integer = length(orblocat),
	period = Bool[1, 1, 0],
	rcut::Real, rcutmat::AbstractMatrix{<:Real} = fill(rcut, (norb, norb)),
	αmax::Real = 1.0, αmaxmat::AbstractMatrix{<:Real} = fill(αmax, (norb, norb)),
	r₀::Real = 0, ϵ₁::Real = 1, ϵ₂::Real = 1, r₀atol::Real = 0,
	δ::Real = 1e-8, head::Union{Number, Nothing} = nothing,
) where {T}

	println("[Optimize mirror correction: keldysh_gauss]")

	println("In general, it's only valid for bare Coulomb potential.")

	δ ≤ 1 || error("δ should be less than or equal to 1.")
	@assert numorb(U) == norb "number of orbitals in U should be the same as length of orblocat."

	period = Vec3{Bool}(period)
	@assert count(period) == 2 "Exactly two dimensions should be periodic."
	if !period[1]
		xyindex = SVector{2, Int}(2, 3)
		zindex = 1
	elseif !period[2]
		xyindex = SVector{2, Int}(3, 1)
		zindex = 2
	else
		xyindex = SVector{2, Int}(1, 2)
		zindex = 3
	end

	@assert kgrid.kgrid_size[zindex] == 1 "The kgrid should be 2D for keldysh_gauss mirror correction."
	@assert all(iszero, lattice[xyindex, zindex]) "Wrong lattice: the non-periodic direction should be perpendicular to the periodic plane."
	@assert all(iszero, lattice[zindex, xyindex]) "Wrong lattice: the non-periodic direction should be perpendicular to the periodic plane."

	(aimr_frac, aimr_frac_xy, aimr_norm2, aimr_norm_xy, aimU, nU_out_of_rcut, remoteU_range) = _getU_out_of_rcut_2D(U, lattice, orblocat, xyindex, rcutmat)

	αmat, αzmat, residmat = _MirrorCorrection_keldysh_gauss_optim(lattice, kgrid, xyindex, zindex, aimr_frac_xy, aimr_norm_xy, aimU, αmaxmat, r₀, ϵ₁, ϵ₂, r₀atol, δ, head)
	# αmat, αzmat, r₀mat, residmat = _MirrorCorrection_keldysh_gauss_optim_r₀(lattice, kgrid, xyindex, zindex, aimr_frac_xy, aimr_norm_xy, aimU, αmaxmat, ϵ₁, ϵ₂, r₀atol, δ, head)

	println("[Optimize mirror correction: keldysh_gauss end]")

	result = MirrorCorrection(Val(:keldysh_gauss), U, kgrid, lattice, orblocat; period, αmat, αzmat, r₀, ϵ₁, ϵ₂, r₀atol, δ, head)

	return (Uwithlr = result.Uwithlr, Umirror = result.Umirror, Ushort = result.Ushort,
		αmat = αmat, αzmat = αzmat, residmat = residmat, remoteU_range = remoteU_range, nU_out_of_rcut = nU_out_of_rcut)
end
function MirrorCorrection(::Val{:keldysh_gauss}, U::AbstractReciprocalHoppings{T}, kgrid::RedKgrid,
	lattice::Lattice, orblocat::Vector{<:ReducedCoordinates}; norb::Integer = length(orblocat),
	period = Bool[1, 1, 0],
	α::Real = 1.0, αmat::AbstractMatrix{<:Real} = fill(α, (norb, norb)),
	αz::Real = 1.0, αzmat::AbstractMatrix{<:Real} = fill(αz, (norb, norb)),
	r₀::Real = 0, ϵ₁::Real = 1, ϵ₂::Real = 1, r₀atol::Real = 0,
	δ::Real = 1e-8, head::Union{Number, Nothing} = nothing,
) where {T}

	println("[Mirror correction: keldysh_gauss]")

	δ ≤ 1 || error("δ should be less than or equal to 1.")
	@assert numorb(U) == norb "number of orbitals in U should be the same as length of orblocat."

	period = Vec3{Bool}(period)
	if !period[1]
		xyindex = SVector{2, Int}(2, 3)
		zindex = 1
	elseif !period[2]
		xyindex = SVector{2, Int}(3, 1)
		zindex = 2
	else
		xyindex = SVector{2, Int}(1, 2)
		zindex = 3
	end

	@assert kgrid.kgrid_size[zindex] == 1 "The kgrid should be 2D for keldysh_gauss mirror correction."
	@assert all(iszero, lattice[xyindex, zindex]) "Wrong lattice: the non-periodic direction should be perpendicular to the periodic plane."
	@assert all(iszero, lattice[zindex, xyindex]) "Wrong lattice: the non-periodic direction should be perpendicular to the periodic plane."

	nk = length(kgrid)
	rlattice = reciprocal(lattice)

	lattice2D = SMatrix{2, 2, Float64, 4}(lattice[xyindex, xyindex])
	rlattice2D = SMatrix{2, 2, Float64, 4}(rlattice[xyindex, xyindex])
	kgridsize2D = SVector{2, Int}(kgrid.kgrid_size[xyindex])

	𝐚 = lattice2D[:, 1]
	𝐛 = lattice2D[:, 2]
	S = abs(𝐚[1] * 𝐛[2] - 𝐚[2] * 𝐛[1])
	halfRz = lattice[zindex, zindex] / 2

	# Create Ggrid for mirror correction.
	allkG, allkG_norm2 = _PeriodCoulombKeldyshGauss_kG(rlattice, kgrid, xyindex, maximum(αmat), δ)
	nkG = length(allkG)
	println("number of k+G: ", nkG * 2 + 1, "\nnumber of positive k+G: ", nkG)

	io_lock = ReentrantLock()
	LinearAlgebra.BLAS.set_num_threads(1)

	Umirror = deepcopy(U)
	Ushort = deepcopy(U)
	tasks = Vector{Task}(undef, norb * (norb + 1) ÷ 2)
	n = 0
	for j in 2:norb, i in 1:(j-1)
		n += 1
		tasks[n] = Threads.@spawn begin
			φkG, φr = _MirrorCorrection_keldysh_gauss_φcor(lattice2D, rlattice2D, kgridsize2D, nk,
				allkG, allkG_norm2, nkG, αmat[i, j], αzmat[i, j], r₀, ϵ₁, ϵ₂, S, halfRz, r₀atol, head)
			#(i, j)
			dorb = orblocat[j] - orblocat[i]
			U_ij = hops(Umirror)[i, j]
			Ushort_ij = hops(Ushort)[i, j]
			for idx in eachindex(U_ij)
				hop = U_ij[idx]
				r_frac = path(hop) + dorb
				r_frac = r_frac[xyindex]
				φkG_r = φkG(r_frac)
				φr_r = φr(r_frac)
				U_ij[idx] = similar(hop, value(hop) - φkG_r + φr_r)
				Ushort_ij[idx] = similar(hop, value(hop) - φkG_r)
			end
			#(j, i)
			dorb = orblocat[i] - orblocat[j]
			U_ji = hops(Umirror)[j, i]
			Ushort_ji = hops(Ushort)[j, i]
			for idx in eachindex(U_ji)
				hop = U_ji[idx]
				r_frac = path(hop) + dorb
				r_frac = r_frac[xyindex]
				φkG_r = φkG(r_frac)
				φr_r = φr(r_frac)
				U_ji[idx] = similar(hop, value(hop) - φkG_r + φr_r)
				Ushort_ji[idx] = similar(hop, value(hop) - φkG_r)
			end
			lock(io_lock) do
				@printf("mirror correction for pair (%d, %d): nU = %d\n", i, j, length(U_ij))
				@printf("mirror correction for pair (%d, %d): nU = %d\n", j, i, length(U_ji))
				flush(stdout)
			end
		end
	end
	for i in 1:norb
		n += 1
		tasks[n] = Threads.@spawn begin
			φkG, φr = _MirrorCorrection_keldysh_gauss_φcor(lattice2D, rlattice2D, kgridsize2D, nk,
				allkG, allkG_norm2, nkG, αmat[i, i], αzmat[i, i], r₀, ϵ₁, ϵ₂, S, halfRz, r₀atol, head)
			#(i, i)
			U_ii = hops(Umirror)[i, i]
			Ushort_ii = hops(Ushort)[i, i]
			for idx in eachindex(U_ii)
				hop = U_ii[idx]
				r_frac = path(hop)
				r_frac = r_frac[xyindex]
				φkG_r = φkG(r_frac)
				φr_r = φr(r_frac)
				U_ii[idx] = similar(hop, value(hop) - φkG_r + φr_r)
				Ushort_ii[idx] = similar(hop, value(hop) - φkG_r)
			end
			lock(io_lock) do
				@printf("mirror correction for pair (%d, %d): nU = %d\n", i, i, length(U_ii))
				flush(stdout)
			end
		end
	end
	wait.(tasks)
	LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())

	Uwithlr = DirectInterAction(:keldysh_gauss, Umirror, lattice, orblocat; αmat, αzmat, r₀, ϵ₁, ϵ₂, r₀atol, δ)
	@printf("Uwithlr: nG = %d\n", Uwithlr.LR.nG)

	println("[Mirror correction: keldysh_gauss end]")

	return (Uwithlr = Uwithlr, Umirror = Umirror, Ushort = Ushort)
end
function _MirrorCorrection_keldysh_gauss_φcor(lattice2D, rlattice2D, kgridsize2D, nk, allkG, allkG_norm2, nkG, α, αz, r₀, ϵ₁, ϵ₂, S, halfRz, r₀atol, head)
	invnk_time_2 = 2 / nk
	φr = RealCoulomb(Val(:keldysh_gauss); r₀, ϵ₁, ϵ₂, α, αz, atol = r₀atol)
	φK = ReciprocalCoulomb(Val(:keldysh_gauss); r₀, ϵ₁, ϵ₂, α, αz, S, atol = r₀atol)
	φkG_data = map(kG2 -> φK(Val(:k²), kG2; atol = 0.0) * invnk_time_2, allkG_norm2)
	if isnothing(head)
		head_value = φK(Val(:head), rlattice2D, kgridsize2D) / nk
	else
		head_value = head + CoulombScale * (2π / S) * (halfRz / 2 - 1 / (sqrt(π) * α_z))
	end
	φkG(r) = sum(ikG -> φkG_data[ikG] * cospi(2 * (allkG[ikG] ⋅ r)), Base.OneTo(nkG)) + head_value
	φr_frac(r) = φr(lattice2D * r)
	return φkG, φr_frac
end
function _MirrorCorrection_keldysh_gauss_optim(lattice, kgrid, xyindex, zindex, aimr_frac_xy, aimr_norm_xy, aimU, αmax, r₀, ϵ₁, ϵ₂, r₀atol, δ, head)

	nk = length(kgrid)
	norb = size(aimU, 1)
	nU = length.(aimU)

	invnk_time_2 = 2 / nk

	rlattice = reciprocal(lattice)

	lattice2D = SMatrix{2, 2, Float64, 4}(lattice[xyindex, xyindex])
	rlattice2D = SMatrix{2, 2, Float64, 4}(rlattice[xyindex, xyindex])
	kgridsize2D = SVector{2, Int}(kgrid.kgrid_size[xyindex])

	𝐚 = lattice2D[:, 1]
	𝐛 = lattice2D[:, 2]
	S = abs(𝐚[1] * 𝐛[2] - 𝐚[2] * 𝐛[1])
	halfRz = lattice[zindex, zindex] / 2

	# Create k+G for optimization.
	allkG, allkG_norm2 = _PeriodCoulombKeldyshGauss_kG(rlattice, kgrid, xyindex, maximum(αmax), δ)
	nkG = length(allkG)
	println("number of k+G: ", nkG * 2 + 1, "\nnumber of positive k+G: ", nkG)

	αmat = Matrix{Float64}(undef, norb, norb)
	αzmat = Matrix{Float64}(undef, norb, norb)
	residmat = Matrix{Float64}(undef, norb, norb)
	isconvergence = falses(norb, norb)

	io_lock = ReentrantLock()

	pair = Vector{Tuple{Int, Int}}(undef, 0)
	sizehint!(pair, norb * (norb + 1) ÷ 2)
	for j in 1:norb, i in 1:j
		push!(pair, (i, j))
	end

	# 保留首次循环，优化区间从0开始。
	for (i, j) in pair
		aimU_ij = aimU[i, j]
		aimr_frac_ij = aimr_frac_xy[i, j]
		aimr_norm2_ij = aimr_norm_xy[i, j]
		αmax_ij = αmax[i, j]
		nU_ij = nU[i, j]

		# shared array for a pair.
		φkG = Vector{Float64}(undef, nkG)
		phase = Matrix{Float64}(undef, nkG, nU_ij) # the bottleneck.
		start_time = time()
		LinearAlgebra.BLAS.set_num_threads(1)
		Threads.@threads for iU in Base.OneTo(nU_ij)
			r = aimr_frac_ij[iU]
			for (ikG, kG) in enumerate(allkG)
				@inbounds phase[ikG, iU] = cospi(2 * (kG ⋅ r)) * invnk_time_2
			end
		end
		LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())
		phasetime = time() - start_time

		result = _MirrorCorrection_keldysh_gauss_optim_ij(rlattice2D, kgridsize2D, nk, allkG_norm2, φkG, phase,
			aimr_norm2_ij, aimU_ij, 0.0, αmax_ij, 0.1, 0.1,
			r₀, ϵ₁, ϵ₂, S, halfRz, r₀atol, head)

		αmat[i, j], αzmat[i, j] = result.minimizer[1], result.minimizer[2]
		residmat[i, j] = result.minimum
		αmat[j, i], αzmat[j, i] = αmat[i, j], αzmat[i, j]
		i ≠ j && (residmat[j, i] = 0.0)
		isconvergence[i, j] = !isapprox(αmat[i, j], αmax_ij; rtol = 1e-2, atol = 1e-3)
		isconvergence[j, i] = isconvergence[i, j]
		lock(io_lock) do
			println("phase size: ($nkG, $nU_ij)")
			println("phase time: $(round(phasetime, digits=4))s")
			@printf("optimize α for pair (%d, %d): α = %.4f, αz = %.4f, resid = %e, steps = %d, runtime = %.2f s, convergence = %s\n",
				i, j, αmat[i, j], αzmat[i, j], residmat[i, j], result.iterations, result.time_run, isconvergence[i, j])
			flush(stdout)
		end
	end

	while !all(isconvergence)
		@printf("unconverged pairs:\n")
		for (i, j) in pair
			if !isconvergence[i, j]
				@printf("(%d, %d): α = %14.8f, αz = %14.8f \n", i, j, αmat[i, j], αzmat[i, j])
			end
		end
		@printf("increase αmax and optimize again.\n")

		# increase the updound of α.
		αmax .*= √2
		# Create k+G for optimization.
		allkG, allkG_norm2 = _PeriodCoulombKeldyshGauss_kG(rlattice, kgrid, xyindex, maximum(αmax), δ)
		nkG = length(allkG)
		println("number of k+G: ", nkG * 2 + 1, "\nnumber of positive k+G: ", nkG)

		# shared array.
		φkG = Vector{Float64}(undef, nkG)
		nUmax = maximum(nU[.!isconvergence])
		phase = Matrix{Float64}(undef, nkG, nUmax)
		flush(stdout)

		for (i, j) in pair
			if isconvergence[i, j]
				continue
			end

			aimU_ij = aimU[i, j]
			aimr_frac_ij = aimr_frac_xy[i, j]
			aimr_norm2_ij = aimr_norm_xy[i, j]
			αmax_ij = αmax[i, j]
			nU_ij = nU[i, j]

			α_ij, αz_ij = αmat[i, j], αzmat[i, j]

			start_time = time()
			LinearAlgebra.BLAS.set_num_threads(1)
			Threads.@threads for iU in Base.OneTo(nU_ij)
				r = aimr_frac_ij[iU]
				for (ikG, kG) in enumerate(allkG)
					@inbounds phase[ikG, iU] = cospi(2 * (kG ⋅ r)) * invnk_time_2
				end
			end
			LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())
			phasetime = time() - start_time
			println("phase size: ($nkG, $(nU_ij))")
			println("phase time: $(round(phasetime, digits=4))s")

			result = _MirrorCorrection_keldysh_gauss_optim_ij(rlattice2D, kgridsize2D, nk, allkG_norm2, φkG, phase,
				aimr_norm2_ij, aimU_ij, αmax_ij/2, αmax_ij, α_ij, αz_ij,
				r₀, ϵ₁, ϵ₂, S, halfRz, r₀atol, head)

			αmat[i, j], αzmat[i, j] = result.minimizer[1], result.minimizer[2]
			residmat[i, j] = result.minimum
			αmat[j, i], αzmat[j, i] = αmat[i, j], αzmat[i, j]
			i ≠ j && (residmat[j, i] = 0.0)
			isconvergence[i, j] = !isapprox(αmat[i, j], αmax_ij; rtol = 1e-2, atol = 1e-3)
			isconvergence[j, i] = isconvergence[i, j]
			@printf("optimize α for pair (%d, %d): α = %.4f, αz = %.4f, resid = %e, steps = %d, runtime = %.2f s, convergence = %s\n",
				i, j, αmat[i, j], αzmat[i, j], residmat[i, j], result.iterations, result.time_run, isconvergence[i, j])
			flush(stdout)
		end
	end

	return αmat, αzmat, residmat
end
function _MirrorCorrection_keldysh_gauss_optim_ij(rlattice2D, kgridsize2D, nk, allkG_norm2, φkG, phase,
	aimr_norm_ij, aimU_ij, αmin_ij, αmax_ij, α0_ij, αz0_ij,
	r₀, ϵ₁, ϵ₂, S, halfRz, r₀atol, head)

	function loss(α)
		φK = ReciprocalCoulomb(Val(:keldysh_gauss); r₀, ϵ₁, ϵ₂, α = α[1], αz = α[2], S, atol = r₀atol)
		map(enumerate(allkG_norm2)) do (ikG, kG2)
			φkG[ikG] = φK(Val(:k²), kG2; atol = 0.0)
		end

		if isnothing(head)
			head_value = φK(Val(:head), rlattice2D, kgridsize2D) / nk
		else
			head_value = head + CoulombScale * (2π / S) * (halfRz / 2 - 1 / (sqrt(π) * α_z))
		end

		resid = sum(enumerate(aimU_ij)) do (iU, v)
			abs2(v - φkG ⋅ view(phase, :, iU) - head_value)
			# abs2(v - φkG ⋅ view(phase, :, iU) - head_value) / aimr_norm_ij[iU]
		end
		return resid
	end

	α0 = [α0_ij, αz0_ij]
	df = TwiceDifferentiable(loss, α0)

	lα = [αmin_ij, 0.0]
	uα = [αmax_ij, Inf]
	dfc = TwiceDifferentiableConstraints(lα, uα)

	return optimize(df, dfc, α0, IPNewton())
end
function _MirrorCorrection_keldysh_gauss_optim_r₀(lattice, kgrid, xyindex, zindex, aimr_frac_xy, aimr_norm_xy, aimU, αmax, ϵ₁, ϵ₂, r₀atol, δ, head)

	nk = length(kgrid)
	norb = size(aimU, 1)
	nU = length.(aimU)

	invnk_time_2 = 2 / nk

	rlattice = reciprocal(lattice)

	lattice2D = SMatrix{2, 2, Float64, 4}(lattice[xyindex, xyindex])
	rlattice2D = SMatrix{2, 2, Float64, 4}(rlattice[xyindex, xyindex])
	kgridsize2D = SVector{2, Int}(kgrid.kgrid_size[xyindex])

	𝐚 = lattice2D[:, 1]
	𝐛 = lattice2D[:, 2]
	S = abs(𝐚[1] * 𝐛[2] - 𝐚[2] * 𝐛[1])
	halfRz = lattice[zindex, zindex] / 2

	# Create k+G for optimization.
	allkG, allkG_norm2 = _PeriodCoulombKeldyshGauss_kG(rlattice, kgrid, xyindex, maximum(αmax), δ)
	nkG = length(allkG)
	println("number of k+G: ", nkG * 2 + 1, "\nnumber of positive k+G: ", nkG)

	αmat = Matrix{Float64}(undef, norb, norb)
	αzmat = Matrix{Float64}(undef, norb, norb)
	r₀mat = Matrix{Float64}(undef, norb, norb)
	residmat = Matrix{Float64}(undef, norb, norb)
	isconvergence = falses(norb, norb)

	io_lock = ReentrantLock()

	pair = Vector{Tuple{Int, Int}}(undef, 0)
	sizehint!(pair, norb * (norb + 1) ÷ 2)
	for j in 1:norb, i in 1:j
		push!(pair, (i, j))
	end

	# 保留首次循环，优化区间从0开始。
	for (i, j) in pair
		aimU_ij = aimU[i, j]
		aimr_frac_ij = aimr_frac_xy[i, j]
		aimr_norm2_ij = aimr_norm_xy[i, j]
		αmax_ij = αmax[i, j]
		nU_ij = nU[i, j]

		# shared array for a pair.
		φkG = Vector{Float64}(undef, nkG)
		phase = Matrix{Float64}(undef, nkG, nU_ij) # the bottleneck.
		start_time = time()
		LinearAlgebra.BLAS.set_num_threads(1)
		Threads.@threads for iU in Base.OneTo(nU_ij)
			r = aimr_frac_ij[iU]
			for (ikG, kG) in enumerate(allkG)
				@inbounds phase[ikG, iU] = cospi(2 * (kG ⋅ r)) * invnk_time_2
			end
		end
		LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())
		phasetime = time() - start_time

		result = _MirrorCorrection_keldysh_gauss_optim_r₀_ij(rlattice2D, kgridsize2D, nk, allkG_norm2, φkG, phase,
			aimr_norm2_ij, aimU_ij, 0.0, αmax_ij, 0.1, 0.1, 1.0,
			ϵ₁, ϵ₂, S, halfRz, r₀atol, head)

		αmat[i, j], αzmat[i, j], r₀mat[i, j] = result.minimizer[1], result.minimizer[2], result.minimizer[3]
		residmat[i, j] = result.minimum
		αmat[j, i], αzmat[j, i], r₀mat[j, i] = αmat[i, j], αzmat[i, j], r₀mat[i, j]
		i ≠ j && (residmat[j, i] = 0.0)
		isconvergence[i, j] = !isapprox(αmat[i, j], αmax_ij; rtol = 1e-2, atol = 1e-3)
		isconvergence[j, i] = isconvergence[i, j]
		lock(io_lock) do
			println("phase size: ($nkG, $nU_ij)")
			println("phase time: $(round(phasetime, digits=4))s")
			@printf("optimize α for pair (%d, %d): α = %.4f, αz = %.4f, r₀ = %.4f, resid = %e, steps = %d, runtime = %.2f s, convergence = %s\n",
				i, j, αmat[i, j], αzmat[i, j], r₀mat[i, j], residmat[i, j], result.iterations, result.time_run, isconvergence[i, j])
			flush(stdout)
		end
	end

	while !all(isconvergence)
		@printf("unconverged pairs:\n")
		for (i, j) in pair
			if !isconvergence[i, j]
				@printf("(%d, %d): α = %14.8f, αz = %14.8f, r₀ = %14.8f \n", i, j, αmat[i, j], αzmat[i, j], r₀mat[i, j])
			end
		end
		@printf("increase αmax and optimize again.\n")

		# increase the updound of α.
		αmax .*= √2
		# Create k+G for optimization.
		allkG, allkG_norm2 = _PeriodCoulombKeldyshGauss_kG(rlattice, kgrid, xyindex, maximum(αmax), δ)
		nkG = length(allkG)
		println("number of k+G: ", nkG * 2 + 1, "\nnumber of positive k+G: ", nkG)

		# shared array.
		φkG = Vector{Float64}(undef, nkG)
		nUmax = maximum(nU[.!isconvergence])
		phase = Matrix{Float64}(undef, nkG, nUmax)
		flush(stdout)

		for (i, j) in pair
			if isconvergence[i, j]
				continue
			end

			aimU_ij = aimU[i, j]
			aimr_frac_ij = aimr_frac_xy[i, j]
			aimr_norm2_ij = aimr_norm_xy[i, j]
			αmax_ij = αmax[i, j]
			nU_ij = nU[i, j]

			α_ij, αz_ij, r₀_ij = αmat[i, j], αzmat[i, j], r₀mat[i, j]

			start_time = time()
			LinearAlgebra.BLAS.set_num_threads(1)
			Threads.@threads for iU in Base.OneTo(nU_ij)
				r = aimr_frac_ij[iU]
				for (ikG, kG) in enumerate(allkG)
					@inbounds phase[ikG, iU] = cospi(2 * (kG ⋅ r)) * invnk_time_2
				end
			end
			LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())
			phasetime = time() - start_time
			println("phase size: ($nkG, $(nU_ij))")
			println("phase time: $(round(phasetime, digits=4))s")

			result = _MirrorCorrection_keldysh_gauss_optim_r₀_ij(rlattice2D, kgridsize2D, nk, allkG_norm2, φkG, phase,
				aimr_norm2_ij, aimU_ij, αmax_ij/2, αmax_ij, α_ij, αz_ij, r₀_ij,
				ϵ₁, ϵ₂, S, halfRz, r₀atol, head)

			αmat[i, j], αzmat[i, j], r₀mat[i, j] = result.minimizer[1], result.minimizer[2], result.minimizer[3]
			residmat[i, j] = result.minimum
			αmat[j, i], αzmat[j, i], r₀mat[j, i] = αmat[i, j], αzmat[i, j], r₀mat[i, j]
			i ≠ j && (residmat[j, i] = 0.0)
			isconvergence[i, j] = !isapprox(αmat[i, j], αmax_ij; rtol = 1e-2, atol = 1e-3)
			isconvergence[j, i] = isconvergence[i, j]
			@printf("optimize α for pair (%d, %d): α = %.4f, αz = %.4f, r₀ = %.4f, resid = %e, steps = %d, runtime = %.2f s, convergence = %s\n",
				i, j, αmat[i, j], αzmat[i, j], r₀mat[i, j], residmat[i, j], result.iterations, result.time_run, isconvergence[i, j])
			flush(stdout)
		end
	end

	return αmat, αzmat, r₀mat, residmat
end
function _MirrorCorrection_keldysh_gauss_optim_r₀_ij(rlattice2D, kgridsize2D, nk, allkG_norm2, φkG, phase,
	aimr_norm_ij, aimU_ij, αmin_ij, αmax_ij, α0_ij, αz0_ij, r₀0,
	ϵ₁, ϵ₂, S, halfRz, r₀atol, head)

	function loss(αr₀)
		φK = ReciprocalCoulomb(Val(:keldysh_gauss); r₀ = αr₀[3], ϵ₁, ϵ₂, α = αr₀[1], αz = αr₀[2], S, atol = r₀atol)
		map(enumerate(allkG_norm2)) do (ikG, kG2)
			φkG[ikG] = φK(Val(:k²), kG2; atol = 0.0)
		end

		if isnothing(head)
			head_value = φK(Val(:head), rlattice2D, kgridsize2D) / nk
		else
			head_value = head + CoulombScale * (2π / S) * (halfRz / 2 - 1 / (sqrt(π) * α_z))
		end

		resid = sum(enumerate(aimU_ij)) do (iU, v)
			abs2(v - φkG ⋅ view(phase, :, iU) - head_value)
			# abs2(v - φkG ⋅ view(phase, :, iU) - head_value) / aimr_norm_ij[iU]
		end
		return resid
	end

	αr₀0 = [α0_ij, αz0_ij, r₀0]
	df = TwiceDifferentiable(loss, αr₀0)

	lαr₀ = [αmin_ij, 0.0, 0.0]
	uαr₀ = [αmax_ij, Inf, Inf]
	dfc = TwiceDifferentiableConstraints(lαr₀, uαr₀)

	return optimize(df, dfc, αr₀0, IPNewton())
end
