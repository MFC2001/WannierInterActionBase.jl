function Base.read(epsmat::AbstractString, ::Type{epsmat_bgw};
	vcoul = joinpath(dirname(epsmat), "vcoul"),
	eps0mat = joinpath(dirname(epsmat), "eps0mat.h5"),
	ωindex::Union{Integer, AbstractVector{<:Integer}} = 1,
	eps0_ωindex = :default,
)
	vcoul = read(vcoul, vcoul_bgw)
	epsmat = h5open(epsmat, "r")

	# qpoints of epsmat
	qirre = map(ReducedCoordinates{Float64}, eachcol(epsmat["eps_header/qpoints/qpts"][]))
	all(isone, epsmat["eps_header/qpoints/qpt_done"][]) || error("Some epsmat isn't be calculated.")
	# qgrid
	qgrid = RedKgrid(MonkhorstPack(epsmat["eps_header/qpoints/qgrid"][]))
	# fractionalize qirre
	qirre = kpoints_Float_to_Rational(qirre, qgrid; atol = 1e-6)

	# ωgrid
	if ωindex isa Number
		ωindex = [ωindex]
	end
	nω = length(ωindex)
	ωgrid = epsmat["/eps_header/freqs/freqs"][]
	ωgrid = complex.(ωgrid[1, :], ωgrid[2, :])
	ωgrid = ωgrid[ωindex]

	# Glist
	Glist_qirre = map(ReducedCoordinates{Int}, eachcol(epsmat["mf_header/gspace/components"][]))
	nG_qirre = Int.(epsmat["eps_header/gspace/nmtx"][])
	Gidx_qirre = epsmat["eps_header/gspace/gind_eps2rho"][]
	Glist_qirre = [Glist_qirre[Int.(Gidx_qirre[1:nG_iq, iq])] for (iq, nG_iq) in enumerate(nG_qirre)]

	# epsmat
	mat = epsmat["mats/matrix"][]
	mat = mat[:, :, :, :, 1, :]
	matrix_flavor = epsmat["eps_header/flavor"][]
	if matrix_flavor == 1
		mat = mat[1, :, :, :, :]
	elseif matrix_flavor == 2
		mat = complex.(mat[1, :, :, :, :], mat[2, :, :, :, :])
	else
		error("Wrong matrix_flavor: $matrix_flavor.")
	end

	nqirre = length(qirre)
	matrices = Matrix{Matrix{eltype(mat)}}(undef, nqirre, nω)
	# note: 排除掉单位矩阵部分。
	for iω in Base.OneTo(nω), iq in 1:nqirre
		matrices[iq, iω] = mat[1:nG_qirre[iq], 1:nG_qirre[iq], ωindex[iω], iq] - I
	end

	#crystal
	alat = epsmat["mf_header/crystal/alat"][] * bohr_radius
	lattice = Lattice(epsmat["mf_header/crystal/avec"][] * alat)
	rlattice = reciprocal(lattice) # ReciprocalLattice(epsmat["mf_header/crystal/bvec"][] * bohr_radius)
	# ldot = parent(lattice)
	# ldot = Mat3(transpose(ldot) * ldot) # Mat3(epsmat["mf_header/crystal/adot"][] * bohr_radius^2)
	rldot = parent(rlattice)
	rldot = Mat3(transpose(rldot) * rldot) # Mat3(epsmat["mf_header/crystal/bdot"][] * bohr_radius^2)
	# Ω = det(parent(lattice)) # epsmat["mf_header/crystal/celvol"][] * bohr_radius^3
	# rΩ = det(parent(rlattice)) # epsmat["mf_header/crystal/recvol"][] * bohr_radius^3

	#symmetry
	nsym = epsmat["mf_header/symmetry/ntran"][]
	symmat = epsmat["mf_header/symmetry/mtrx"][]
	symtran = epsmat["mf_header/symmetry/tnp"][]
	syms = Vector{SymOp}(undef, nsym)
	for isym in Base.OneTo(nsym)
		syms[isym] = SymOp(transpose(symmat[:, :, isym]), symtran[:, isym])
	end

	close(epsmat)

	# symmetrize epsmat
	Threads.@threads for iq in Base.OneTo(nqirre)
		q = qirre[iq]
		iq_v = findfirst(q_v->norm(q_v - q) ≤ 1e-7, vcoul.qpoints)
		_symmetrize_epsmat_bgw!(view(matrices, iq, :), vcoul.vals[iq_v])
	end

	# eps0mat
	if eps0mat isa AbstractString
		eps0mat = [eps0mat]
	end
	neps0mat = length(eps0mat)
	if eps0_ωindex == :default
		eps0_ωindex = fill!(ωindex, neps0mat)
	elseif any(x->length(x) ≠ nω, eps0_ωindex)
		error("Mismatched length of ωindex of epsmat and eps0mat!")
	end

	# matrices_q0, q0, Glist_q0
	matrices_q0 = Matrix{Matrix{eltype(mat)}}(undef, neps0mat, nω)
	q0 = Vector{CartesianCoordinates{Float64}}(undef, neps0mat)
	Glist_q0 = Vector{Vector{ReducedCoordinates{Int}}}(undef, neps0mat)
	for ieps0mat in 1:neps0mat
		eps0mat = h5open(eps0mat[ieps0mat], "r")
		q0_frac = ReducedCoordinates{Float64}(eps0mat["eps_header/qpoints/qpts"][])
		q0[ieps0mat] = rlattice * q0_frac
		Glist_q0_ieps0mat = map(ReducedCoordinates{Int}, eachcol(eps0mat["mf_header/gspace/components"][]))
		nG_q0 = eps0mat["eps_header/gspace/nmtx"][]
		nG_q0 = Int(nG_q0[1])
		Gidx_q0 = eps0mat["eps_header/gspace/gind_eps2rho"][]
		Glist_q0[ieps0mat] = Glist_q0_ieps0mat[Int.(Gidx_q0[1:nG_q0, 1])]
		# epsmat
		mat = eps0mat["mats/matrix"][]
		mat = mat[:, :, :, :, 1, :]
		matrix_flavor = eps0mat["eps_header/flavor"][]
		if matrix_flavor == 1
			mat = mat[1, :, :, :, :]
		elseif matrix_flavor == 2
			mat = complex.(mat[1, :, :, :, :], mat[2, :, :, :, :])
		else
			error("Wrong matrix_flavor: $matrix_flavor.")
		end
		# note: 排除掉单位矩阵部分。
		for iω in Base.OneTo(nω)
			matrices_q0[ieps0mat, iω] = mat[1:nG_q0, 1:nG_q0, eps0_ωindex[ieps0mat][iω], 1] - I
		end
		close(eps0mat)
		# symmetrize epsmat
		iq_v = findfirst(q_v->all(≤(1e-7), q_v - q0_frac), vcoul.qpoints)
		_symmetrize_epsmat_bgw!(view(matrices_q0, ieps0mat, :), vcoul.vals[iq_v])
	end

	# use a same Glist for all q0.
	Glist_q0_unique = intersect(Glist_q0...)
	Glist_q0_unique = sort!(Glist_q0_unique; by = G->G ⋅ (rldot * G))
	for ieps0mat in 1:neps0mat
		G2idx = Dict(G => iG for (iG, G) in enumerate(Glist_q0[ieps0mat]))
		Gidx = map(G -> G2idx[G], Glist_q0_unique)
		for iω in Base.OneTo(nω)
			matrices_q0[ieps0mat, iω] = matrices_q0[ieps0mat, iω][Gidx, Gidx]
		end
	end

	return DielectricPWBasis(matrices, qirre, Glist_qirre, matrices_q0, q0, Glist_q0_unique, qgrid, ωgrid, syms, lattice)
end
function _symmetrize_epsmat_bgw!(mat_qω_inv_iq, vcoul_iq)
	vcoul_sqrt_iq = sqrt.(vcoul_iq)
	vcoul_mat = [v_iG2 / v_iG1 for v_iG1 in vcoul_sqrt_iq, v_iG2 in vcoul_sqrt_iq]
	for mat in mat_qω_inv_iq
		mat .*= vcoul_mat
	end
	return mat_qω_inv_iq
end
