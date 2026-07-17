# 使用3D高斯分布的库伦势进行镜像修正，适用于3D系统。
# 暂不支持多极修正
function MirrorCorrection(::Val{:gauss_optim}, U::AbstractReciprocalHoppings{T}, kgrid::RedKgrid,
	lattice::Lattice, orblocat::Vector{<:ReducedCoordinates}; norb::Integer = length(orblocat),
	rcut::Real, rcutmat::AbstractMatrix{<:Real} = fill(rcut, (norb, norb)),
	rcutmin::Real = 8, rcutminmat::AbstractMatrix{<:Real} = fill(rcutmin, (norb, norb)),
	αmax::Union{Real, Nothing} = nothing,
	αmaxmat::AbstractMatrix{<:Real} = isnothing(αmax) ? 4.5 ./ rcutminmat : fill(αmax, (norb, norb)),
	ϵ::Real = 1, δ::Real = 1e-8, head::Union{Number, Nothing} = nothing,
) where {T}

	println("[Optimize mirror correction: gauss]")

	println("In general, it's only valid for bare Coulomb potential.")

	δ ≤ 1 || error("δ should be less than or equal to 1.")
	@assert numorb(U) == norb "number of orbitals in U should be the same as length of orblocat."

	(aimr_frac, aimr_norm2, aimU, nU_out_of_rcut, remoteU_range) = _getU_out_of_rcut_3D(U, lattice, orblocat, rcutmat)

	αmat, residmat = _MirrorCorrection_gauss_optim(lattice, kgrid, aimr_frac, aimr_norm2, aimU, αmaxmat, ϵ, δ, head)

	println("[Optimize mirror correction: gauss end]")

	result = MirrorCorrection(Val(:gauss), U, kgrid, lattice, orblocat; αmat, ϵ, δ, head)

	return (Uwithlr = result.Uwithlr, Umirror = result.Umirror, Ushort = result.Ushort,
		αmat = αmat, residmat = residmat, remoteU_range = remoteU_range, nU_out_of_rcut = nU_out_of_rcut)
end
function MirrorCorrection(::Val{:gauss}, U::AbstractReciprocalHoppings{T}, kgrid::RedKgrid,
	lattice::Lattice, orblocat::Vector{<:ReducedCoordinates}; norb::Integer = length(orblocat),
	rcut::Real = 1.0, rcutmat::AbstractMatrix{<:Real} = fill(rcut, (norb, norb)),
	α::Union{Real, Nothing} = nothing,
	αmat::AbstractMatrix{<:Real} = isnothing(α) ? 4.5 ./ rcutmat : fill(α, (norb, norb)),
	ϵ::Real = 1, δ::Real = 1e-8, head::Union{Number, Nothing} = nothing,
) where {T}

	println("[Mirror correction: gauss]")

	δ ≤ 1 || error("δ should be less than or equal to 1.")
	@assert numorb(U) == norb "number of orbitals in U should be the same as length of orblocat."

	nk = length(kgrid)

	Ω = abs(det(parent(lattice)))
	rlattice = reciprocal(lattice)

	# Create Ggrid for mirror correction.
	allkG, allkG_norm2 = _PeriodCoulombGauss_kG(rlattice, kgrid, maximum(αmat), δ)
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
			φkG, φr = _MirrorCorrection_gauss_φcor(lattice, rlattice, kgrid.kgrid_size, nk,
				allkG, allkG_norm2, nkG, αmat[i, j], ϵ, Ω, head)
			#(i, j)
			dorb = orblocat[j] - orblocat[i]
			U_ij = hops(Umirror)[i, j]
			Ushort_ij = hops(Ushort)[i, j]
			for idx in eachindex(U_ij)
				hop = U_ij[idx]
				r_frac = path(hop) + dorb
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
			φkG, φr = _MirrorCorrection_gauss_φcor(lattice, rlattice, kgrid.kgrid_size, nk,
				allkG, allkG_norm2, nkG, αmat[i, i], ϵ, Ω, head)
			#(i, i)
			U_ii = hops(Umirror)[i, i]
			Ushort_ii = hops(Ushort)[i, i]
			for idx in eachindex(U_ii)
				hop = U_ii[idx]
				r_frac = path(hop)
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

	Uwithlr = DirectInterAction(:gauss, Umirror, lattice, orblocat; αmat, δ, ϵ)
	@printf("Uwithlr: nG = %d\n", Uwithlr.LR.nG)

	println("[Mirror correction: gauss end]")

	return (Uwithlr = Uwithlr, Umirror = Umirror, Ushort = Ushort)
end
function _MirrorCorrection_gauss_φcor(lattice, rlattice, kgridsize, nk, allkG, allkG_norm2, nkG, α, ϵ, Ω, head)
	invnk_time_2 = 2 / nk
	φr = RealCoulomb(Val(:gauss); ϵ, α)
	φK = ReciprocalCoulomb(Val(:gauss); ϵ, α, Ω)
	φkG_data = map(kG2 -> φK(Val(:k²), kG2; atol = 0.0) * invnk_time_2, allkG_norm2)
	if isnothing(head)
		head_value = φK(Val(:head), rlattice, kgridsize) / nk
	else
		head_value = head - CoulombScale * (4π / (ϵ * Ω)) / (4 * α * α)
	end
	φkG(r) = sum(ikG -> φkG_data[ikG] * cospi(2 * (allkG[ikG] ⋅ r)), Base.OneTo(nkG)) + head_value
	φr_frac(r) = φr(lattice * r)
	return φkG, φr_frac
end
function _MirrorCorrection_gauss_optim(lattice, kgrid, aimr_frac, aimr_norm2, aimU, αmax, ϵ, δ, head)

	nk = length(kgrid)
	norb = size(aimU, 1)
	nU = length.(aimU)

	invnk_time_2 = 2 / nk

	Ω = abs(det(parent(lattice)))
	rlattice = reciprocal(lattice)

	# Create k+G for optimization.
	allkG, allkG_norm2 = _PeriodCoulombGauss_kG(rlattice, kgrid, maximum(αmax), δ)
	nkG = length(allkG)
	println("number of k+G: ", nkG * 2 + 1, "\nnumber of positive k+G: ", nkG)

	αmat = Matrix{Float64}(undef, norb, norb)
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
		aimr_frac_ij = aimr_frac[i, j]
		aimr_norm2_ij = aimr_norm2[i, j]
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

		result = _MirrorCorrection_gauss_optim_ij(rlattice, kgrid.kgrid_size, nk,
			allkG_norm2, φkG, phase, aimr_norm2_ij, aimU_ij, 0.0, αmax_ij, ϵ, Ω, head)

		αmat[i, j] = result.minimizer
		residmat[i, j] = result.minimum
		αmat[j, i] = αmat[i, j]
		i ≠ j && (residmat[j, i] = 0.0)
		isconvergence[i, j] = !isapprox(αmat[i, j], αmax_ij; rtol = 1e-2, atol = 1e-3)
		isconvergence[j, i] = isconvergence[i, j]
		lock(io_lock) do
			println("phase size: ($nkG, $nU_ij)")
			println("phase time: $(round(phasetime, digits=4))s")
			@printf("optimize α for pair (%d, %d): α = %.4f, resid = %e, steps = %d, runtime = %.2f s, convergence = %s\n",
				i, j, αmat[i, j], residmat[i, j], result.iterations, result.time_run, isconvergence[i, j])
			flush(stdout)
		end
	end

	while !all(isconvergence)
		@printf("unconverged pairs:\n")
		for (i, j) in pair
			if !isconvergence[i, j]
				@printf("(%d, %d): α = %14.8f \n", i, j, αmat[i, j])
			end
		end
		@printf("increase αmax and optimize again.\n")

		# increase the updound of α.
		αmax .*= √2
		# Create k+G for optimization.
		allkG, allkG_norm2 = _PeriodCoulombGauss_kG(rlattice, kgrid, maximum(αmax), δ)
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
			aimr_frac_ij = aimr_frac[i, j]
			aimr_norm2_ij = aimr_norm2[i, j]
			αmax_ij = αmax[i, j]
			nU_ij = nU[i, j]

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

			result = _MirrorCorrection_gauss_optim_ij(rlattice, kgrid.kgrid_size, nk,
				allkG_norm2, φkG, phase, aimr_norm2_ij, aimU_ij, αmax_ij/2, αmax_ij, ϵ, Ω, head)

			αmat[i, j] = result.minimizer
			residmat[i, j] = result.minimum
			αmat[j, i] = αmat[i, j]
			i ≠ j && (residmat[j, i] = 0.0)
			isconvergence[i, j] = !isapprox(αmat[i, j], αmax_ij; rtol = 1e-2, atol = 1e-3)
			isconvergence[j, i] = isconvergence[i, j]
			@printf("optimize α for pair (%d, %d): α = %.4f, resid = %e, steps = %d, runtime = %.2f s, convergence = %s\n",
				i, j, αmat[i, j], residmat[i, j], result.iterations, result.time_run, isconvergence[i, j])
			flush(stdout)
		end
	end

	return αmat, residmat
end
function _MirrorCorrection_gauss_optim_ij(rlattice, kgridsize, nk, allkG_norm2, φkG, phase,
	aimr_norm2_ij, aimU_ij, αmin_ij, αmax_ij, ϵ, Ω, head)

	function loss(α)
		φK = ReciprocalCoulomb(Val(:gauss); ϵ, α, Ω)
		map(enumerate(allkG_norm2)) do (ikG, kG2)
			φkG[ikG] = φK(Val(:k²), kG2; atol = 0.0)
		end

		if isnothing(head)
			head_value = φK(Val(:head), rlattice, kgridsize) / nk
		else
			head_value = head - CoulombScale * (4π / (ϵ * Ω)) / (4 * α * α)
		end

		# α增大，此求和形式的ϕk：近距离增大，远距离减小。
		resid = sum(enumerate(aimU_ij)) do (iU, v)
			abs2(v - φkG ⋅ view(phase, :, iU) - head_value)
			# abs2(v - φkG ⋅ view(phase, :, iU) - head_value) / aimr_norm2_ij[iU]
		end
		return resid
	end

	return optimize(loss, αmin_ij, αmax_ij, Brent())
end
