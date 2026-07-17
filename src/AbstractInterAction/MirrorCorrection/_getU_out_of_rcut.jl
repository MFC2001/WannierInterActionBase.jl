function _getU_out_of_rcut_3D(U::AbstractReciprocalHoppings{T}, lattice, orblocat, rcutmat) where {T}
	ldot = parent(lattice)
	ldot = Mat3(transpose(ldot) * ldot)

	norb = numorb(U)
	aimr_frac = Matrix{Vector{ReducedCoordinates{Float64}}}(undef, norb, norb)
	aimr_norm2 = Matrix{Vector{Float64}}(undef, norb, norb)
	aimU = Matrix{Vector{T}}(undef, norb, norb)
	nU_out_of_rcut = Matrix{Int}(undef, norb, norb)
	remoteU_range = Matrix{Tuple{Float64, Float64}}(undef, norb, norb)

	for j in 1:norb, i in 1:norb
		dorb = orblocat[j] - orblocat[i]
		aimr_frac_ij = Vector{ReducedCoordinates{Float64}}(undef, 0)
		aimr_norm2_ij = Vector{Float64}(undef, 0)
		aimU_ij = Vector{T}(undef, 0)
		rcut2 = rcutmat[i, j]^2
		for hop in hops(U)[i, j]
			r_frac = path(hop) + dorb
			r2 = r_frac ⋅ (ldot * r_frac)
			if r2 > rcut2
				push!(aimr_frac_ij, r_frac)
				push!(aimr_norm2_ij, r2)
				push!(aimU_ij, value(hop))
			end
		end
		if isempty(aimU_ij)
			error("no interaction between orbital $i and $j within rcut $(rcutmat[i, j]).")
		end
		nU_out_of_rcut[i, j] = length(aimU_ij)
		aimr_frac[i, j] = aimr_frac_ij
		aimr_norm2[i, j] = aimr_norm2_ij
		aimU[i, j] = aimU_ij
		(r2min, r2max) = extrema(aimr_norm2_ij)
		rmin, rmax = sqrt(r2min), sqrt(r2max)
		remoteU_range[i, j] = (rmin, rmax)
		@printf("pair (%d, %d) have %d U in range [%.2f, %.2f] out of rcut %.4f.\n",
			i, j, nU_out_of_rcut[i, j], rmin, rmax, rcutmat[i, j])
	end
	flush(stdout)
	return aimr_frac, aimr_norm2, aimU, nU_out_of_rcut, remoteU_range
end
function _getU_out_of_rcut_2D(U::AbstractReciprocalHoppings{T}, lattice, orblocat, xyindex, rcutmat) where {T}
	ldot = parent(lattice)
	ldot = Mat3(transpose(ldot) * ldot)
	lattice2D = SMatrix{2, 2, Float64, 4}(lattice[xyindex, xyindex])
	ldot_xy = SMatrix{2, 2, Float64, 4}(transpose(lattice2D) * lattice2D)

	norb = numorb(U)
	aimr_frac = Matrix{Vector{ReducedCoordinates{Float64}}}(undef, norb, norb)
	aimr_frac_xy = Matrix{Vector{SVector{2, Float64}}}(undef, norb, norb)
	aimr_norm2 = Matrix{Vector{Float64}}(undef, norb, norb)
	aimr_norm_xy = Matrix{Vector{Float64}}(undef, norb, norb)
	aimU = Matrix{Vector{T}}(undef, norb, norb)
	nU_out_of_rcut = Matrix{Int}(undef, norb, norb)
	remoteU_range = Matrix{Tuple{Float64, Float64}}(undef, norb, norb)

	for j in 1:norb, i in 1:norb
		dorb = orblocat[j] - orblocat[i]
		aimr_frac_ij = Vector{ReducedCoordinates{Float64}}(undef, 0)
		aimr_frac_xy_ij = Vector{SVector{2, Float64}}(undef, 0)
		aimr_norm2_ij = Vector{Float64}(undef, 0)
		aimr_norm_xy_ij = Vector{Float64}(undef, 0)
		aimU_ij = Vector{T}(undef, 0)
		rcut2 = rcutmat[i, j]^2
		for hop in hops(U)[i, j]
			r_frac = path(hop) + dorb
			r_frac_xy = r_frac[xyindex]
			r2_xy = r_frac_xy ⋅ (ldot_xy * r_frac_xy)
			if r2_xy > rcut2
				push!(aimr_frac_ij, r_frac)
				push!(aimr_frac_xy_ij, r_frac_xy)
				push!(aimr_norm2_ij, r_frac ⋅ (ldot * r_frac))
				push!(aimr_norm_xy_ij, sqrt(r2_xy))
				push!(aimU_ij, value(hop))
			end
		end
		if isempty(aimU_ij)
			error("no interaction between orbital $i and $j within rcut $(rcutmat[i, j]).")
		end
		nU_out_of_rcut[i, j] = length(aimU_ij)
		aimr_frac[i, j] = aimr_frac_ij
		aimr_frac_xy[i, j] = aimr_frac_xy_ij
		aimr_norm2[i, j] = aimr_norm2_ij
		aimr_norm_xy[i, j] = aimr_norm_xy_ij
		aimU[i, j] = aimU_ij
		(r2min, r2max) = extrema(aimr_norm2_ij)
		rmin, rmax = sqrt(r2min), sqrt(r2max)
		remoteU_range[i, j] = (rmin, rmax)
		@printf("pair (%d, %d) have %d U in range [%.2f, %.2f] out of rcut %.4f.\n",
			i, j, nU_out_of_rcut[i, j], rmin, rmax, rcutmat[i, j])
	end
	flush(stdout)
	return aimr_frac, aimr_frac_xy, aimr_norm2, aimr_norm_xy, aimU, nU_out_of_rcut, remoteU_range
end
