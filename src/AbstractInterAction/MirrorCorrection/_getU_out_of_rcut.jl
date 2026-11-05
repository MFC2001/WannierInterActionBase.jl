function _getU_out_of_rcut_3D(U::HRFull{T}, lattice, orblocat, rcutmat) where {T}
    ldot = parent(lattice)
    ldot = Mat3(transpose(ldot) * ldot)
    rcut2mat = rcutmat .* rcutmat

    norb = numorb(U)
    dorb = [orblocat[j] - orblocat[i] for i in 1:norb, j in 1:norb]
    aimr_frac = [Vector{ReducedCoordinates{Float64}}(undef, 0) for _ in 1:norb, _ in 1:norb]
    aimr_norm2 = [Vector{Float64}(undef, 0) for _ in 1:norb, _ in 1:norb]
    aimU = [Vector{T}(undef, 0) for _ in 1:norb, _ in 1:norb]
    aimU_ndegen = [Vector{Int}(undef, 0) for _ in 1:norb, _ in 1:norb]
    nU_out_of_rcut = Matrix{Int}(undef, norb, norb)
    remoteU_range = Matrix{Tuple{Float64,Float64}}(undef, norb, norb)

    for (i, j, iR, val) in enumerate(U)
        r_frac = U.R[iR] + dorb[i, j]
        r2 = r_frac ⋅ (ldot * r_frac)
        if r2 > rcut2mat[i, j]
            push!(aimr_frac[i, j], r_frac)
            push!(aimr_norm2[i, j], r2)
            push!(aimU[i, j], val)
            push!(aimU_ndegen[i, j], U.ndegen[iR])
        end
    end
    for j in 1:norb, i in 1:norb
        if isempty(aimU[i, j])
            error("no interaction between orbital $i and $j within rcut $(rcutmat[i, j]).")
        end
        nU_out_of_rcut[i, j] = length(aimU[i, j])
        (r2min, r2max) = extrema(aimr_norm2[i, j])
        rmin, rmax = sqrt(r2min), sqrt(r2max)
        remoteU_range[i, j] = (rmin, rmax)
        @printf("pair (%d, %d) have %d U in range [%.2f, %.2f] out of rcut %.4f.\n",
            i, j, nU_out_of_rcut[i, j], rmin, rmax, rcutmat[i, j])
    end
    flush(stdout)
    return aimr_frac, aimr_norm2, aimU, aimU_ndegen, nU_out_of_rcut, remoteU_range
end
function _getU_out_of_rcut_2D(U::HRFull{T}, lattice, orblocat, xyindex, rcutmat) where {T}
    ldot = parent(lattice)
    ldot = Mat3(transpose(ldot) * ldot)
    lattice2D = SMatrix{2,2,Float64,4}(lattice[xyindex, xyindex])
    ldot_xy = SMatrix{2,2,Float64,4}(transpose(lattice2D) * lattice2D)
    rcut2mat = rcutmat .* rcutmat

    norb = numorb(U)
    dorb = [orblocat[j] - orblocat[i] for i in 1:norb, j in 1:norb]
    aimr_frac = [Vector{ReducedCoordinates{Float64}}(undef, 0) for _ in 1:norb, _ in 1:norb]
    aimr_frac_xy = [Vector{SVector{2,Float64}}(undef, 0) for _ in 1:norb, _ in 1:norb]
    aimr_norm2 = [Vector{Float64}(undef, 0) for _ in 1:norb, _ in 1:norb]
    aimr_norm_xy = [Vector{Float64}(undef, 0) for _ in 1:norb, _ in 1:norb]
    aimU = [Vector{T}(undef, 0) for _ in 1:norb, _ in 1:norb]
    aimU_ndegen = [Vector{Int}(undef, 0) for _ in 1:norb, _ in 1:norb]
    nU_out_of_rcut = Matrix{Int}(undef, norb, norb)
    remoteU_range = Matrix{Tuple{Float64,Float64}}(undef, norb, norb)

    for (i, j, iR, val) in enumerate(U)
        r_frac = U.R[iR] + dorb[i, j]
        r_frac_xy = r_frac[xyindex]
        r2_xy = r_frac_xy ⋅ (ldot_xy * r_frac_xy)
        if r2_xy > rcut2mat[i, j]
            push!(aimr_frac[i, j], r_frac)
            push!(aimr_frac_xy[i, j], r_frac_xy)
            push!(aimr_norm2[i, j], r_frac ⋅ (ldot * r_frac))
            push!(aimr_norm_xy[i, j], sqrt(r2_xy))
            push!(aimU[i, j], val)
            push!(aimU_ndegen[i, j], U.ndegen[iR])
        end
    end
    for j in 1:norb, i in 1:norb
        if isempty(aimU[i, j])
            error("no interaction between orbital $i and $j within rcut $(rcutmat[i, j]).")
        end
        nU_out_of_rcut[i, j] = length(aimU[i, j])
        (rmin, rmax) = extrema(aimr_norm_xy[i, j])
        remoteU_range[i, j] = (rmin, rmax)
        @printf("pair (%d, %d) have %d U in range [%.2f, %.2f] out of rcut %.4f.\n",
            i, j, nU_out_of_rcut[i, j], rmin, rmax, rcutmat[i, j])
    end
    flush(stdout)
    return aimr_frac, aimr_frac_xy, aimr_norm2, aimr_norm_xy, aimU, aimU_ndegen, nU_out_of_rcut, remoteU_range
end
