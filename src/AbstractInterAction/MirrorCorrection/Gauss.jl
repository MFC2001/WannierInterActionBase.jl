# 使用3D高斯分布的库伦势进行镜像修正，适用于3D系统。
# 暂不支持多极修正
function MirrorCorrection(::Val{:gauss_optim}, U::HRFull{T}, kgrid::RedKgrid,
    lattice::Lattice, orblocat::Vector{<:ReducedCoordinates}; norb::Integer=length(orblocat),
    rcut::Real, rcutmat::AbstractMatrix{<:Real}=fill(rcut, (norb, norb)),
    rcutmin::Real=8, rcutminmat::AbstractMatrix{<:Real}=fill(rcutmin, (norb, norb)),
    αmax::Union{Real,Nothing}=nothing,
    αmaxmat::AbstractMatrix{<:Real}=isnothing(αmax) ? 4.5 ./ rcutminmat : fill(αmax, (norb, norb)),
    ϵ::Real=1, ϵmat::AbstractMatrix{<:Real}=fill(ϵ, (norb, norb)), ϵ_optim::Bool=false,
    δ::Real=1e-8, head::Union{Number,Nothing}=nothing,
) where {T}

    println("[Optimize mirror correction: gauss]")

    println("In general, it's only valid for bare Coulomb potential.")
    println("Because the dielectric constant always has small deviations.")

    δ ≤ 1 || error("δ should be less than or equal to 1.")
    @assert numorb(U) == norb "number of orbitals in U should be the same as length of orblocat."

    (aimr_frac, aimr_norm2, aimU, aimU_ndegen, nU_out_of_rcut, remoteU_range) = _getU_out_of_rcut_3D(U, lattice, orblocat, rcutmat)

    αmat, ϵmat, residmat = _MirrorCorrection_gauss_optim(lattice, kgrid, aimr_frac, aimr_norm2, aimU, aimU_ndegen, αmaxmat, ϵmat, ϵ_optim, δ, head)

    println("[Optimize mirror correction: gauss end]")

    result = MirrorCorrection(Val(:gauss), U, kgrid, lattice, orblocat; αmat, ϵmat, δ, head)

    return (Uwithlr=result.Uwithlr, Umirror=result.Umirror, Ushort=result.Ushort,
        αmat=αmat, ϵmat=ϵmat, residmat=residmat, remoteU_range=remoteU_range, nU_out_of_rcut=nU_out_of_rcut)
end
function MirrorCorrection(::Val{:gauss}, U::HRFull{T}, kgrid::RedKgrid,
    lattice::Lattice, orblocat::Vector{<:ReducedCoordinates}; norb::Integer=length(orblocat),
    rcut::Real=1.0, rcutmat::AbstractMatrix{<:Real}=fill(rcut, (norb, norb)),
    α::Union{Real,Nothing}=nothing,
    αmat::AbstractMatrix{<:Real}=isnothing(α) ? 4.5 ./ rcutmat : fill(α, (norb, norb)),
    ϵ::Real=1, ϵmat::AbstractMatrix{<:Real}=fill(ϵ, (norb, norb)),
    δ::Real=1e-8, head::Union{Number,Nothing}=nothing,
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
    println("number of k+G: ", nkG * 2 + 1)
    println("number of positive k+G: ", nkG)
    println("number of U: ", length(U))

    dorb = [orblocat[j] - orblocat[i] for i in 1:norb, j in 1:norb]
    φkG = Matrix{Function}(undef, norb, norb)
    φr = Matrix{Function}(undef, norb, norb)
    for j in 1:norb, i in 1:norb
        φkG[i, j], φr[i, j] = _MirrorCorrection_gauss_φcor(lattice, rlattice, kgrid.kgrid_size, nk,
            allkG, allkG_norm2, nkG, αmat[i, j], ϵmat[i, j], Ω, head)
    end

    Umirror = deepcopy(U)
    Ushort = deepcopy(U)
    nthreads_blas = LinearAlgebra.BLAS.get_num_threads()
    LinearAlgebra.BLAS.set_num_threads(1)
    Threads.@threads for ih in Base.OneTo(length(U))
        i, j, iR = index(U, ih)
        val = U[ih]
        r_frac = U.R[iR] + dorb[i, j]
        φkG_r = φkG[i, j](r_frac)
        φr_r = φr[i, j](r_frac)
        Umirror[ih] = val - φkG_r + φr_r
        Ushort[ih] = val - φkG_r
    end
    LinearAlgebra.BLAS.set_num_threads(nthreads_blas)

    SR_min = fill(Inf, norb, norb)
    for (i, j, _, v) in enumerate(Ushort)
        v_abs = abs(v)
        if v_abs < SR_min[i, j]
            SR_min[i, j] = v_abs
        end
    end
    for j in Base.OneTo(norb), i in Base.OneTo(norb)
        @printf("(%d, %d): minimum of short-range part is %.6f\n", i, j, SR_min[i, j])
    end

    Uwithlr = DirectInterAction(:gauss, Umirror, lattice, orblocat; αmat, ϵmat, δ)
    @printf("nG of long-range part is %d\n", Uwithlr.LR.nG)

    println("[Mirror correction: gauss end]")
    flush(stdout)

    return (Uwithlr=Uwithlr, Umirror=Umirror, Ushort=Ushort)
end
function _MirrorCorrection_gauss_φcor(lattice, rlattice, kgridsize, nk, allkG, allkG_norm2, nkG, α, ϵ, Ω, head)
    invnk_times_2 = 2 / nk
    φr = RealCoulomb(Val(:gauss); ϵ, α)
    φK = ReciprocalCoulomb(Val(:gauss); ϵ, α, Ω)
    φkG_data = map(kG2 -> φK(Val(:k²), kG2; atol=0.0) * invnk_times_2, allkG_norm2)
    if isnothing(head)
        head_value = φK(Val(:head), rlattice, kgridsize) / nk
    else
        head_value = head - CoulombScale * (4π / (ϵ * Ω)) / (4 * α * α)
    end
    φkG(r) = sum(ikG -> φkG_data[ikG] * cospi(2 * (allkG[ikG] ⋅ r)), Base.OneTo(nkG)) + head_value
    φr_frac(r) = φr(lattice * r)
    return φkG, φr_frac
end
function _MirrorCorrection_gauss_optim(lattice, kgrid, aimr_frac, aimr_norm2, aimU, aimU_ndegen, αmax, ϵmat, ϵ_optim, δ, head)

    nk = length(kgrid)
    norb = size(aimU, 1)
    nU = length.(aimU)

    invnk_times_2 = 2 / nk

    Ω = abs(det(parent(lattice)))
    rlattice = reciprocal(lattice)

    # Create k+G for optimization.
    allkG, allkG_norm2 = _PeriodCoulombGauss_kG(rlattice, kgrid, maximum(αmax), δ)
    nkG = length(allkG)
    println("number of k+G: ", nkG * 2 + 1, "\nnumber of positive k+G: ", nkG)

    αmin = zeros(Float64, norb, norb)
    αmat = Matrix{Float64}(undef, norb, norb)
    ϵmat = Float64.(ϵmat)
    residmat = Matrix{Float64}(undef, norb, norb)
    isconvergence = falses(norb, norb)

    pair = Vector{Tuple{Int,Int}}(undef, 0)
    sizehint!(pair, norb * (norb + 1) ÷ 2)
    for j in 1:norb, i in 1:j
        push!(pair, (i, j))
    end

    # 保留首次循环，优化区间从0开始。
    for (i, j) in pair
        aimr_frac_ij = aimr_frac[i, j]
        nU_ij = nU[i, j]

        # shared array for a pair.
        φkG = Vector{Float64}(undef, nkG)
        phase = Matrix{Float64}(undef, nkG, nU_ij) # the bottleneck.
        start_time = time()
        LinearAlgebra.BLAS.set_num_threads(1)
        # Threads.@threads for iU in Base.OneTo(nU_ij)
        #     r = aimr_frac_ij[iU]
        #     for (ikG, kG) in enumerate(allkG)
        #         @inbounds phase[ikG, iU] = cospi(2 * (kG ⋅ r)) * invnk_times_2
        #     end
        # end
        Threads.@threads for I in CartesianIndices(phase)
            (ikG, iU) = Tuple(I)
            @inbounds phase[ikG, iU] = cospi(2 * (allkG[ikG] ⋅ aimr_frac_ij[iU])) * invnk_times_2
        end
        LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())
        phasetime = time() - start_time
        println("phase size: ($nkG, $nU_ij)")
        println("phase time: $(round(phasetime, digits=4))s")

        result = _MirrorCorrection_gauss_optim_ij(Val(ϵ_optim), i, j, rlattice, kgrid.kgrid_size, nk,
            allkG_norm2, φkG, phase, aimr_norm2, aimU, aimU_ndegen, αmin, αmax, 0.1 * αmax[i, j], ϵmat[i, j], Ω, head)

        αmat[i, j] = result.α
        ϵmat[i, j] = result.ϵ
        residmat[i, j] = result.resid
        αmat[j, i] = result.α
        i ≠ j && (residmat[j, i] = 0.0)
        isconvergence[i, j] = result.isconvergence
        isconvergence[j, i] = result.isconvergence
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
        αmin .= αmax ./ 2
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

            aimr_frac_ij = aimr_frac[i, j]
            nU_ij = nU[i, j]

            start_time = time()
            LinearAlgebra.BLAS.set_num_threads(1)
            Threads.@threads for I in CartesianIndices(phase)
                (ikG, iU) = Tuple(I)
                @inbounds phase[ikG, iU] = cospi(2 * (allkG[ikG] ⋅ aimr_frac_ij[iU])) * invnk_times_2
            end
            LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())
            phasetime = time() - start_time
            println("phase size: ($nkG, $(nU_ij))")
            println("phase time: $(round(phasetime, digits=4))s")

            result = _MirrorCorrection_gauss_optim_ij(Val(ϵ_optim), i, j, rlattice, kgrid.kgrid_size, nk,
                allkG_norm2, φkG, phase, aimr_norm2, aimU, aimU_ndegen, αmin, αmax, αmat[i, j], ϵmat[i, j], Ω, head)

            αmat[i, j] = result.α
            ϵmat[i, j] = result.ϵ
            residmat[i, j] = result.resid
            αmat[j, i] = result.α
            i ≠ j && (residmat[j, i] = 0.0)
            isconvergence[i, j] = result.isconvergence
            isconvergence[j, i] = result.isconvergence
        end
    end

    return αmat, ϵmat, residmat
end
function _MirrorCorrection_gauss_optim_ij(::Val{false}, i, j, rlattice, kgridsize, nk, allkG_norm2, φkG, phase,
    aimr_norm2, aimU, aimU_ndegen, αmin, αmax, α0_ij, ϵ0_ij, Ω, head)

    aimr_norm2_ij = aimr_norm2[i, j]
    aimU_ij = aimU[i, j]
    aimU_ndegen_ij = aimU_ndegen[i, j]
    αmin_ij = αmin[i, j]
    αmax_ij = αmax[i, j]

    function loss(α)
        φK = ReciprocalCoulomb(Val(:gauss); ϵ=ϵ0_ij, α, Ω)
        map(enumerate(allkG_norm2)) do (ikG, kG2)
            φkG[ikG] = φK(Val(:k²), kG2; atol=0.0)
        end

        if isnothing(head)
            head_value = φK(Val(:head), rlattice, kgridsize) / nk
        else
            head_value = head - CoulombScale * (4π / (ϵ0_ij * Ω)) / (4 * α * α)
        end

        # α增大，此求和形式的ϕk：近距离增大，远距离减小。
        resid = sum(enumerate(aimU_ij)) do (iU, v)
            abs2(v - φkG ⋅ view(phase, :, iU) - head_value) / aimU_ndegen_ij[iU]
            # abs2(v - φkG ⋅ view(phase, :, iU) - head_value) / aimr_norm2_ij[iU]
        end
        return resid
    end

    α_tol = 1e-10
    αmin_ij = αmin_ij < α_tol ? α_tol : αmin_ij
    result = optimize(loss, αmin_ij, αmax_ij, Brent())

    α = result.minimizer
    resid = result.minimum

    if !Optim.converged(result)
        error_string = @sprintf("optimization not converged at pair (%d, %d): α = %.4f, resid = %e, upper=%.4f, iterations=%d",
            i, j, α, resid, αmax_ij, result.iterations)
        error(error_string)
    end

    isconvergence = !isapprox(α, αmax_ij; rtol=1e-2, atol=1e-3)
    @printf("optimize α for pair (%d, %d): α = %.4f, resid = %e, steps = %d, runtime = %.2f s, convergence = %s\n",
        i, j, α, resid, result.iterations, result.time_run, isconvergence)
    flush(stdout)

    return (α=α, ϵ=ϵ0_ij, resid=resid, isconvergence=isconvergence)
end
function _MirrorCorrection_gauss_optim_ij(::Val{true}, i, j, rlattice, kgridsize, nk, allkG_norm2, φkG, phase,
    aimr_norm2, aimU, aimU_ndegen, αmin, αmax, α0_ij, ϵ0_ij, Ω, head)

    aimr_norm2_ij = aimr_norm2[i, j]
    aimU_ij = aimU[i, j]
    aimU_ndegen_ij = aimU_ndegen[i, j]
    αmin_ij = αmin[i, j]
    αmax_ij = αmax[i, j]

    function loss(αϵ)

        α = exp(αϵ[1]) + αmin_ij
        ϵ = αϵ[2]

        φK = ReciprocalCoulomb(Val(:gauss); ϵ, α, Ω)
        map(enumerate(allkG_norm2)) do (ikG, kG2)
            φkG[ikG] = φK(Val(:k²), kG2; atol=0.0)
        end

        if isnothing(head)
            head_value = φK(Val(:head), rlattice, kgridsize) / nk
        else
            head_value = head - CoulombScale * (4π / (ϵ * Ω)) / (4 * α * α)
        end

        # α增大，此求和形式的ϕk：近距离增大，远距离减小。
        resid = sum(enumerate(aimU_ij)) do (iU, v)
            abs2(v - φkG ⋅ view(phase, :, iU) - head_value) / aimU_ndegen_ij[iU]
            # abs2(v - φkG ⋅ view(phase, :, iU) - head_value) / aimr_norm2_ij[iU]
        end
        return resid
    end

    w = αmax_ij - αmin_ij
    α0_ij = α0_ij < αmin_ij + 0.1 * w ? αmin_ij + 0.1 * w : α0_ij
    ϵ0_ij = ϵ0_ij < 1.01 ? 1.01 : ϵ0_ij
    αϵ0 = [log(α0_ij - αmin_ij), ϵ0_ij]
    df = TwiceDifferentiable(loss, αϵ0)

    lαϵ = [-Inf, 1.0]
    uαϵ = [log(αmax_ij - αmin_ij), Inf]
    dfc = TwiceDifferentiableConstraints(lαϵ, uαϵ)

    result = optimize(df, dfc, αϵ0, IPNewton())

    α = exp(result.minimizer[1]) + αmin_ij
    ϵ = result.minimizer[2]
    resid = result.minimum

    if !Optim.converged(result)
        error_string = @sprintf("optimization not converged at pair (%d, %d): α = %.4f, ϵ = %.4f, resid = %e, upper=%.4f, iterations=%d",
            i, j, α, ϵ, resid, αmax_ij, result.iterations)
        error(error_string)
    end

    isconvergence = !isapprox(α, αmax_ij; rtol=1e-2, atol=1e-3)
    @printf("optimize α for pair (%d, %d): α = %.4f, ϵ = %.4f, resid = %e, steps = %d, runtime = %.2f s, convergence = %s\n",
        i, j, α, ϵ, resid, result.iterations, result.time_run, isconvergence)
    flush(stdout)

    return (α=α, ϵ=ϵ, resid=resid, isconvergence=isconvergence)
end
