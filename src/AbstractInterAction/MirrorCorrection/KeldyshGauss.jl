# 适用于2D系统。
function MirrorCorrection(::Val{:keldysh_gauss_optim}, U::HRFull{T}, kgrid::RedKgrid,
    lattice::Lattice, orblocat::Vector{<:ReducedCoordinates}; norb::Integer=length(orblocat),
    period=Bool[1, 1, 0],
    rcut::Real, rcutmat::AbstractMatrix{<:Real}=fill(rcut, (norb, norb)),
    αmax::Real=1.0, αmaxmat::AbstractMatrix{<:Real}=fill(αmax, (norb, norb)),
    r₀::Real=0.0, r₀mat::AbstractMatrix{<:Real}=fill(r₀, (norb, norb)), r₀_optim::Bool=false,
    ϵ₁::Real=1, ϵ₂::Real=1, δ::Real=1e-8, head::Union{Number,Nothing}=nothing,
) where {T}

    println("[Optimize mirror correction: keldysh_gauss]")

    δ ≤ 1 || error("δ should be less than or equal to 1.")
    @assert numorb(U) == norb "number of orbitals in U should be the same as length of orblocat."

    period = Vec3{Bool}(period)
    @assert count(period) == 2 "Exactly two dimensions should be periodic."
    if !period[1]
        xyindex = SVector{2,Int}(2, 3)
        zindex = 1
    elseif !period[2]
        xyindex = SVector{2,Int}(3, 1)
        zindex = 2
    else
        xyindex = SVector{2,Int}(1, 2)
        zindex = 3
    end

    @assert kgrid.kgrid_size[zindex] == 1 "The kgrid should be 2D for keldysh_gauss mirror correction."
    @assert all(iszero, lattice[xyindex, zindex]) "Wrong lattice: the non-periodic direction should be perpendicular to the periodic plane."
    @assert all(iszero, lattice[zindex, xyindex]) "Wrong lattice: the non-periodic direction should be perpendicular to the periodic plane."

    (aimr_frac, aimr_frac_xy, aimr_norm2, aimr_norm_xy, aimU, aimU_ndegen, nU_out_of_rcut, remoteU_range) = _getU_out_of_rcut_2D(U, lattice, orblocat, xyindex, rcutmat)

    αmat, αzmat, r₀mat, residmat = _MirrorCorrection_keldysh_gauss_optim(lattice, kgrid, xyindex, zindex, aimr_frac_xy, aimr_norm_xy, aimU, aimU_ndegen, αmaxmat, r₀mat, r₀_optim, ϵ₁, ϵ₂, δ, head)

    println("[Optimize mirror correction: keldysh_gauss end]")

    result = MirrorCorrection(Val(:keldysh_gauss), U, kgrid, lattice, orblocat; period, αmat, αzmat, r₀mat, ϵ₁, ϵ₂, δ, head)

    return (Uwithlr=result.Uwithlr, Umirror=result.Umirror, Ushort=result.Ushort,
        αmat=αmat, αzmat=αzmat, r₀mat=r₀mat, residmat=residmat, remoteU_range=remoteU_range, nU_out_of_rcut=nU_out_of_rcut)
end
function MirrorCorrection(::Val{:keldysh_gauss}, U::HRFull{T}, kgrid::RedKgrid,
    lattice::Lattice, orblocat::Vector{<:ReducedCoordinates}; norb::Integer=length(orblocat),
    period=Bool[1, 1, 0],
    α::Real=1.0, αmat::AbstractMatrix{<:Real}=fill(α, (norb, norb)),
    αz::Real=1.0, αzmat::AbstractMatrix{<:Real}=fill(αz, (norb, norb)),
    r₀::Real=0, r₀mat::AbstractMatrix{<:Real}=fill(r₀, (norb, norb)),
    ϵ₁::Real=1, ϵ₂::Real=1, δ::Real=1e-8, head::Union{Number,Nothing}=nothing,
) where {T}

    println("[Mirror correction: keldysh_gauss]")

    δ ≤ 1 || error("δ should be less than or equal to 1.")
    @assert numorb(U) == norb "number of orbitals in U should be the same as length of orblocat."

    period = Vec3{Bool}(period)
    if !period[1]
        xyindex = SVector{2,Int}(2, 3)
        zindex = 1
    elseif !period[2]
        xyindex = SVector{2,Int}(3, 1)
        zindex = 2
    else
        xyindex = SVector{2,Int}(1, 2)
        zindex = 3
    end

    @assert kgrid.kgrid_size[zindex] == 1 "The kgrid should be 2D for keldysh_gauss mirror correction."
    @assert all(iszero, lattice[xyindex, zindex]) "Wrong lattice: the non-periodic direction should be perpendicular to the periodic plane."
    @assert all(iszero, lattice[zindex, xyindex]) "Wrong lattice: the non-periodic direction should be perpendicular to the periodic plane."

    nk = length(kgrid)
    rlattice = reciprocal(lattice)

    lattice2D = SMatrix{2,2,Float64,4}(lattice[xyindex, xyindex])
    rlattice2D = SMatrix{2,2,Float64,4}(rlattice[xyindex, xyindex])
    kgridsize2D = SVector{2,Int}(kgrid.kgrid_size[xyindex])

    𝐚 = lattice2D[:, 1]
    𝐛 = lattice2D[:, 2]
    S = abs(𝐚[1] * 𝐛[2] - 𝐚[2] * 𝐛[1])
    halfRz = lattice[zindex, zindex] / 2

    # Create Ggrid for mirror correction.
    allkG, allkG_norm2 = _PeriodCoulombKeldyshGauss_kG(rlattice, kgrid, xyindex, maximum(αmat), δ)
    nkG = length(allkG)
    println("number of k+G: ", nkG * 2 + 1, "\nnumber of positive k+G: ", nkG)
    println("number of U: ", length(U))

    dorb = [orblocat[j] - orblocat[i] for i in 1:norb, j in 1:norb]
    φkG = Matrix{Function}(undef, norb, norb)
    φr = Matrix{Function}(undef, norb, norb)
    for j in 1:norb, i in 1:norb
        φkG[i, j], φr[i, j] = _MirrorCorrection_keldysh_gauss_φcor(lattice2D, rlattice2D, kgridsize2D, nk,
            allkG, allkG_norm2, nkG, αmat[i, j], αzmat[i, j], r₀mat[i, j], ϵ₁, ϵ₂, S, halfRz, head)
    end

    Umirror = deepcopy(U)
    Ushort = deepcopy(U)
    nthreads_blas = LinearAlgebra.BLAS.get_num_threads()
    LinearAlgebra.BLAS.set_num_threads(1)
    Threads.@threads for ih in Base.OneTo(length(U))
        i, j, iR = index(U, ih)
        val = U[ih]
        r_frac = U.R[iR] + dorb[i, j]
        r_frac_xy = r_frac[xyindex]
        φkG_r = φkG[i, j](r_frac_xy)
        φr_r = φr[i, j](r_frac_xy)
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

    Uwithlr = DirectInterAction(:keldysh_gauss, Umirror, lattice, orblocat; αmat, αzmat, r₀mat, ϵ₁, ϵ₂, δ)
    @printf("nG of long-range part is %d\n", Uwithlr.LR.nG)

    println("[Mirror correction: keldysh_gauss end]")
    flush(stdout)

    return (Uwithlr=Uwithlr, Umirror=Umirror, Ushort=Ushort)
end
function _MirrorCorrection_keldysh_gauss_φcor(lattice2D, rlattice2D, kgridsize2D, nk, allkG, allkG_norm2, nkG, α, αz, r₀, ϵ₁, ϵ₂, S, halfRz, head)
    invnk_times_2 = 2 / nk
    φr = RealCoulomb(Val(:keldysh_gauss); r₀, ϵ₁, ϵ₂, α, αz)
    φK = ReciprocalCoulomb(Val(:keldysh_gauss); r₀, ϵ₁, ϵ₂, α, αz, S)
    φkG_data = map(kG2 -> φK(Val(:k²), kG2; atol=0.0) * invnk_times_2, allkG_norm2)
    if isnothing(head)
        head_value = φK(Val(:head), rlattice2D, kgridsize2D) / nk
    else
        head_value = head + CoulombScale * (2π / S) * (halfRz / 2 - 1 / (sqrt(π) * αz))
    end
    φkG(r) = sum(ikG -> φkG_data[ikG] * cospi(2 * (allkG[ikG] ⋅ r)), Base.OneTo(nkG)) + head_value
    φr_frac(r) = φr(lattice2D * r)
    return φkG, φr_frac
end
function _MirrorCorrection_keldysh_gauss_optim(lattice, kgrid, xyindex, zindex, aimr_frac_xy, aimr_norm_xy, aimU, aimU_ndegen, αmax, r₀mat, r₀_optim, ϵ₁, ϵ₂, δ, head)

    nk = length(kgrid)
    norb = size(aimU, 1)
    nU = length.(aimU)

    invnk_times_2 = 2 / nk

    rlattice = reciprocal(lattice)

    lattice2D = SMatrix{2,2,Float64,4}(lattice[xyindex, xyindex])
    rlattice2D = SMatrix{2,2,Float64,4}(rlattice[xyindex, xyindex])
    kgridsize2D = SVector{2,Int}(kgrid.kgrid_size[xyindex])

    𝐚 = lattice2D[:, 1]
    𝐛 = lattice2D[:, 2]
    S = abs(𝐚[1] * 𝐛[2] - 𝐚[2] * 𝐛[1])
    halfRz = lattice[zindex, zindex] / 2

    # Create k+G for optimization.
    allkG, allkG_norm2 = _PeriodCoulombKeldyshGauss_kG(rlattice, kgrid, xyindex, maximum(αmax), δ)
    nkG = length(allkG)
    println("number of k+G: ", nkG * 2 + 1, "\nnumber of positive k+G: ", nkG)

    αmin = zeros(Float64, norb, norb)
    αmat = Matrix{Float64}(undef, norb, norb)
    αzmat = Matrix{Float64}(undef, norb, norb)
    r₀mat = Float64.(r₀mat)
    residmat = Matrix{Float64}(undef, norb, norb)
    isconvergence = falses(norb, norb)

    pair = Vector{Tuple{Int,Int}}(undef, 0)
    sizehint!(pair, norb * (norb + 1) ÷ 2)
    for j in 1:norb, i in 1:j
        push!(pair, (i, j))
    end

    # 保留首次循环，优化区间从0开始。
    for (i, j) in pair
        aimr_frac_ij = aimr_frac_xy[i, j]
        nU_ij = nU[i, j]

        # shared array for a pair.
        φkG = Vector{Float64}(undef, nkG)
        phase = Matrix{Float64}(undef, nkG, nU_ij) # the bottleneck.
        start_time = time()
        LinearAlgebra.BLAS.set_num_threads(1)
        Threads.@threads for I in CartesianIndices(phase)
            (ikG, iU) = Tuple(I)
            @inbounds phase[ikG, iU] = cospi(2 * (allkG[ikG] ⋅ aimr_frac_ij[iU])) * invnk_times_2
        end
        LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())
        phasetime = time() - start_time
        println("phase size: ($nkG, $nU_ij)")
        println("phase time: $(round(phasetime, digits=4))s")

        result = _MirrorCorrection_keldysh_gauss_optim_ij(Val(r₀_optim), i, j,
            rlattice2D, kgridsize2D, nk, allkG_norm2, φkG, phase,
            aimr_norm_xy, aimU, aimU_ndegen, αmin, αmax, 0.1 * αmax[i, j], 0.1, r₀mat[i, j],
            ϵ₁, ϵ₂, S, halfRz, head)

        αmat[i, j], αzmat[i, j], r₀mat[i, j] = result.α, result.αz, result.r₀
        residmat[i, j] = result.resid
        αmat[j, i], αzmat[j, i], r₀mat[j, i] = result.α, result.αz, result.r₀
        i ≠ j && (residmat[j, i] = 0.0)
        isconvergence[i, j] = result.isconvergence
        isconvergence[j, i] = result.isconvergence
    end

    while !all(isconvergence)
        @printf("unconverged pairs:\n")
        for (i, j) in pair
            if !isconvergence[i, j]
                @printf("(%d, %d): α = %14.8f, αz = %14.8f, r₀ =  %14.8f \n", i, j, αmat[i, j], αzmat[i, j], r₀mat[i, j])
            end
        end
        @printf("increase αmax and optimize again.\n")

        # increase the updound of α.
        αmax .*= √2
        αmin .= αmax ./ 2
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

            aimr_frac_ij = aimr_frac_xy[i, j]
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

            result = _MirrorCorrection_keldysh_gauss_optim_ij(Val(r₀_optim), i, j,
                rlattice2D, kgridsize2D, nk, allkG_norm2, φkG, phase,
                aimr_norm_xy, aimU, aimU_ndegen, αmin, αmax, αmat[i, j], αzmat[i, j], r₀mat[i, j],
                ϵ₁, ϵ₂, S, halfRz, head)

            αmat[i, j], αzmat[i, j], r₀mat[i, j] = result.α, result.αz, result.r₀
            residmat[i, j] = result.resid
            αmat[j, i], αzmat[j, i], r₀mat[j, i] = result.α, result.αz, result.r₀
            i ≠ j && (residmat[j, i] = 0.0)
            isconvergence[i, j] = result.isconvergence
            isconvergence[j, i] = result.isconvergence
        end
    end

    return αmat, αzmat, r₀mat, residmat
end
function _MirrorCorrection_keldysh_gauss_optim_ij(::Val{false}, i, j,
    rlattice2D, kgridsize2D, nk, allkG_norm2, φkG, phase,
    aimr_norm, aimU, aimU_ndegen, αmin, αmax, α0_ij, αz0_ij, r₀0_ij,
    ϵ₁, ϵ₂, S, halfRz, head)

    aimr_norm_ij = aimr_norm[i, j]
    aimU_ij = aimU[i, j]
    aimU_ndegen_ij = aimU_ndegen[i, j]
    αmin_ij = αmin[i, j]
    αmax_ij = αmax[i, j]

    function loss(α_in)

        α = exp(α_in[1]) + αmin_ij
        αz = exp(α_in[2])

        φK = ReciprocalCoulomb(Val(:keldysh_gauss); r₀=r₀0_ij, ϵ₁, ϵ₂, α, αz, S)
        map(enumerate(allkG_norm2)) do (ikG, kG2)
            φkG[ikG] = φK(Val(:k²), kG2; atol=0.0)
        end

        if isnothing(head)
            head_value = φK(Val(:head), rlattice2D, kgridsize2D) / nk
        else
            head_value = head + CoulombScale * (2π / S) * (halfRz / 2 - 1 / (sqrt(π) * αz))
        end

        resid = sum(enumerate(aimU_ij)) do (iU, v)
            abs2(v - φkG ⋅ view(phase, :, iU) - head_value) / aimU_ndegen_ij[iU]
            # abs2(v - φkG ⋅ view(phase, :, iU) - head_value) / aimr_norm_ij[iU]
        end
        return resid
    end

    w = αmax_ij - αmin_ij
    α0_ij = α0_ij < αmin_ij + 0.1 * w ? αmin_ij + 0.1 * w : α0_ij
    αz0_ij = αz0_ij < 0.01 ? 0.01 : αz0_ij
    α0 = [log(α0_ij - αmin_ij), log(αz0_ij)]
    df = TwiceDifferentiable(loss, α0)

    lα = [-Inf, -Inf]
    uα = [log(αmax_ij - αmin_ij), Inf]
    dfc = TwiceDifferentiableConstraints(lα, uα)

    result = optimize(df, dfc, α0, IPNewton())

    α = exp(result.minimizer[1]) + αmin_ij
    αz = exp(result.minimizer[2])
    resid = result.minimum

    if !Optim.converged(result)
        error_string = @sprintf("optimization not converged at pair (%d, %d): α = %.4f, αz = %.4f, resid = %e, α_upper=%.4f, iterations=%d",
            i, j, α, αz, resid, αmax_ij, result.iterations)
        error(error_string)
    end

    isconvergence = !isapprox(α, αmax_ij; rtol=1e-2, atol=1e-3)
    @printf("optimize α for pair (%d, %d): α = %.4f, αz = %.4f, resid = %e, steps = %d, runtime = %.2f s, convergence = %s\n",
        i, j, α, αz, resid, result.iterations, result.time_run, isconvergence)
    flush(stdout)

    return (α=α, αz=αz, r₀=r₀0_ij, resid=resid, isconvergence=isconvergence)
end
function _MirrorCorrection_keldysh_gauss_optim_ij(::Val{true}, i, j,
    rlattice2D, kgridsize2D, nk, allkG_norm2, φkG, phase,
    aimr_norm, aimU, aimU_ndegen, αmin, αmax, α0_ij, αz0_ij, r₀0_ij,
    ϵ₁, ϵ₂, S, halfRz, head)

    aimr_norm_ij = aimr_norm[i, j]
    aimU_ij = aimU[i, j]
    aimU_ndegen_ij = aimU_ndegen[i, j]
    αmin_ij = αmin[i, j]
    αmax_ij = αmax[i, j]

    function loss(αr₀)

        α = exp(αr₀[1]) + αmin_ij
        αz = exp(αr₀[2])
        r₀ = exp(αr₀[3])

        φK = ReciprocalCoulomb(Val(:keldysh_gauss); r₀, ϵ₁, ϵ₂, α, αz, S)
        map(enumerate(allkG_norm2)) do (ikG, kG2)
            φkG[ikG] = φK(Val(:k²), kG2; atol=0.0)
        end

        if isnothing(head)
            head_value = φK(Val(:head), rlattice2D, kgridsize2D) / nk
        else
            head_value = head + CoulombScale * (2π / S) * (halfRz / 2 - 1 / (sqrt(π) * αz))
        end

        resid = sum(enumerate(aimU_ij)) do (iU, v)
            abs2(v - φkG ⋅ view(phase, :, iU) - head_value) / aimU_ndegen_ij[iU]
            # abs2(v - φkG ⋅ view(phase, :, iU) - head_value) / aimr_norm_ij[iU]
        end
        return resid
    end

    w = αmax_ij - αmin_ij
    α0_ij = α0_ij < αmin_ij + 0.1 * w ? αmin_ij + 0.1 * w : α0_ij
    αz0_ij = αz0_ij < 0.01 ? 0.01 : αz0_ij
    r₀0_ij = r₀0_ij < 0.01 ? 0.01 : r₀0_ij
    αr₀0 = [log(α0_ij - αmin_ij), log(αz0_ij), log(r₀0_ij)]
    df = TwiceDifferentiable(loss, αr₀0)

    lαr₀ = [-Inf, -Inf, -Inf]
    uαr₀ = [log(αmax_ij - αmin_ij), Inf, Inf]
    dfc = TwiceDifferentiableConstraints(lαr₀, uαr₀)

    result = optimize(df, dfc, αr₀0, IPNewton())

    α = exp(result.minimizer[1]) + αmin_ij
    αz = exp(result.minimizer[2])
    r₀ = exp(result.minimizer[3])
    resid = result.minimum

    if !Optim.converged(result)
        error_string = @sprintf("optimization not converged at pair (%d, %d): α = %.4f, αz = %.4f, r₀ = %.4f, resid = %e, α_upper=%.4f, iterations=%d",
            i, j, α, αz, r₀, resid, αmax_ij, result.iterations)
        error(error_string)
    end

    isconvergence = !isapprox(α, αmax_ij; rtol=1e-2, atol=1e-3)
    @printf("optimize α for pair (%d, %d): α = %.4f, αz = %.4f, r₀ = %.4f, resid = %e, steps = %d, runtime = %.2f s, convergence = %s\n",
        i, j, α, αz, r₀, resid, result.iterations, result.time_run, isconvergence)
    flush(stdout)

    return (α=α, αz=αz, r₀=r₀, resid=resid, isconvergence=isconvergence)
end
