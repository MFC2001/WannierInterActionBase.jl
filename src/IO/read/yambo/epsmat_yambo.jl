function Base.read(epsfolder::AbstractString, ::Type{epsmat_yambo};
	SAVE::AbstractString = joinpath(parentdir(epsfolder), "SAVE"),
	q̂₀ = CartesianCoordinates{Float64}(0, 0, 1),
	format::Symbol = :netcdf4,
	suffix::Union{AbstractString, Nothing} = nothing,
	ωindex::Union{Integer, AbstractVector{<:Integer}} = 1,
	qgrid::Union{MonkhorstPack, Nothing} = nothing, qgrid_atol::Real = 1e-6,
	eps0folder = nothing, eps0_ωindex = :default, eps0_q̂₀ = nothing,
)

	isdir(epsfolder) || error("$epsfolder isn't found!")

	if isnothing(suffix)
		suffixes = ["em1s", "pp", "Xmpa", "em1d"]
		for s in suffixes
			if isfile(joinpath(epsfolder, "ndb." * s))
				suffix = s
				break
			end
		end
	end
	isnothing(suffix) && error("Can't find dielectric data in $epsfolder, please check if them exist or specify the suffix of the file.")
	isdir(SAVE) || error("SAVE folder not found in $(parentdir(epsfolder)).")

	db1 = NCDataset(joinpath(SAVE, "ns.db1"), "r"; format)
	lattice_parameter = db1["LATTICE_PARAMETER"][:]
	rlattice_parameter = 2π ./ lattice_parameter
	lattice = Lattice(transpose(db1["LATTICE_VECTORS"][:, :]))
	rlattice = reciprocal(lattice)
	rldot = transpose(parent(rlattice)) * parent(rlattice)
	# in Å.
	lattice_Å = Lattice(transpose(db1["LATTICE_VECTORS"][:, :] .* bohr_radius))

	symmat = db1["SYMMETRY"][:, :, :]
	nsym = size(symmat, 3)
	syms = Vector{SymOp}(undef, nsym)
	lmat = parent(lattice)
	lmat_inv = inv(lmat)
	for isym in Base.OneTo(nsym)
		S = round.(Int, lmat_inv * symmat[:, :, isym] * lmat)
		syms[isym] = SymOp(S, [0, 0, 0])
	end

	# dimensions = db1["DIMENSIONS"][:]
	# wf_nG = db1["WF_G_COMPONENTS"][:]
	# max_atoms = db1["MAX_ATOMS"][:]
	# number_of_atom_species = db1["number_of_atom_species"][:]
	# exx_fraction = db1["EXX_FRACTION"][:]
	# exx_screening = db1["EXX_SCREENING"][:]
	# mag_syms = db1["mag_syms"][:]
	# GPL_revision = db1["GPL_REVISION"][:]
	# n_atoms = db1["N_ATOMS"][:]
	# atom_pos = db1["ATOM_POS"][:, :, :]
	# atomic_numbers = db1["atomic_numbers"][:]
	# atom_mass = db1["ATOM_MASS"][:]
	# atom_map = db1["ATOM_MAP"][:,:]
	# Glist_rho = db1["G-VECTORS"][:,:]
	# k_irre = db1["K-POINTS"][:,:]
	# eigenvalues = db1["EIGENVALUES"][:,:,:]
	# nG_wave = db1["WFC_NG"][:]
	# Gidx_wave = db1["WFC_GRID"][:,:]
	close(db1)

	eps_head = NCDataset(joinpath(epsfolder, "ndb." * suffix), "r"; format)
	spin_vars = Int.(eps_head["SPIN_VARS"][:])
	if spin_vars[1] == 1 # 应该是自旋通道数
		if spin_vars[2] == 1
			nspin = 1
		elseif spin_vars[2] == 2 # 应该是波函数的旋量分量数
			nspin = 4
		else
			error("Invalid SPIN_VARS: $(spin_vars).")
		end
	elseif spin_vars[1] == 2
		if spin_vars[2] == 1
			nspin = 2
		else
			error("Invalid SPIN_VARS: $(spin_vars).")
		end
	else
		error("Invalid SPIN_VARS: $(spin_vars).")
	end

	nkq = eps_head["HEAD_R_LATT"][:]
	nqirre = Int(nkq[1])
	nq = Int(nkq[2])
	nk_irre = Int(nkq[3])
	nk = Int(nkq[4])

	qirre = eps_head["HEAD_QPT"][:, :]
	qirre = map(Base.OneTo(nqirre)) do iq
		rlattice \ CartesianCoordinates{Float64}(rlattice_parameter .* qirre[iq, :])
	end

	# qgrid
	if qgrid isa MonkhorstPack
		qirre = kpoints_Float_to_Rational(qirre, qgrid)
		qgrid = RedKgrid(qgrid)
	else
		qirre = kpoints_Float_to_Rational(qirre; tol = qgrid_atol)
		allq = Set{ReducedCoordinates{Rational{Int}}}()
		for q in qirre
			for sym in syms
				push!(allq, normalize_kdirect(sym.S * q))
			end
		end
		allq = collect(allq)
		qx = unique(map(q->q[1], allq))
		qy = unique(map(q->q[2], allq))
		qz = unique(map(q->q[3], allq))
		nqx = length(qx)
		nqy = length(qy)
		nqz = length(qz)
		@assert nqx * nqy * nqz == nq "Cannot determine qgrid from q_irre, please check if the eps_data are correct or adjust the `qgrid_atol` or set `qgrid` explicitly."
		qgrid = RedKgrid(MonkhorstPack(nqx, nqy, nqz))
	end

	Glist = eps_head["X_RL_vecs"][:, :]
	nG = size(Glist, 1)
	Glist = map(Base.OneTo(nG)) do iG
		round.(Int, rlattice \ CartesianCoordinates{Float64}(rlattice_parameter .* Glist[iG, :]))
	end

	# head_wf = eps_head["HEAD_WF"][:]
	# cutoff = eps_head["CUTOFF"][:, :]
	# fragmented = eps_head["FRAGMENTED"][:]
	# temperatures = eps_head["TEMPERATURES"][:]
	# gauge = eps_head["GAUGE"][:]
	# x_pars_1 = eps_head["X_PARS_1"][:]
	# x_time_ordering = eps_head["X_Time_ordering"][:]
	# x_tddft_kernel = eps_head["X_TDDFT_KERNEL"][:]
	# x_drued = eps_head["X_DRUDE"][:]
	# x_pars_3 = eps_head["X_PARS_3"][:]
	# x_optical_average = eps_head["X_OPTICAL_AVERAGE"][:]
	close(eps_head)

	# epsmat
	nω = length(ωindex)
	matrices = Matrix{Matrix{ComplexF64}}(undef, nqirre, nω)
	eps = NCDataset(joinpath(epsfolder, "ndb." * suffix * "_fragment_1"), "r"; format)
	# freq_pars_sec_iq = eps["FREQ_PARS_sec_iq1"][:]
	ωgrid = eps["FREQ_sec_iq1"][:, :]
	ωgrid = complex.(ωgrid[1, :], ωgrid[2, :])
	ωgrid = ωgrid[ωindex]
	matrices_iq1 = eps["X_Q_1"][:, :, :, ωindex]
	close(eps)
	for iω in Base.OneTo(nω)
		#note: yambo存储的是不包含单位矩阵部分，且完全对称化的介电函数。
		matrices[1, iω] = complex.(matrices_iq1[1, :, :, iω], matrices_iq1[2, :, :, iω])
	end
	for iq in 2:nqirre
		eps = NCDataset(joinpath(epsfolder, "ndb." * suffix * "_fragment_$(iq)"), "r"; format)
		matrices_iq = eps["X_Q_$(iq)"][:, :, :, ωindex]
		close(eps)
		for iω in Base.OneTo(nω)
			#note: yambo存储的是不包含单位矩阵部分，且完全对称化的介电函数。
			matrices[iq, iω] = complex.(matrices_iq[1, :, :, iω], matrices_iq[2, :, :, iω])
		end
	end

	q0_idx = findfirst(q->iszero(q), qirre)
	q_else = setdiff(eachindex(qirre), q0_idx)
	matrices_qirre = matrices[q_else, :]
	qirre = qirre[q_else]
	Glist_qirre = map(i -> Glist, Base.OneTo(nqirre - 1))

	# other eps for q0
	if isnothing(eps0folder)
		matrices_q0 = matrices[[q0_idx], :]
		q0 = [CartesianCoordinates{Float64}(q̂₀)]
		Glist_q0 = Glist
		return DielectricPWBasis(matrices_qirre, qirre, Glist_qirre, matrices_q0, q0, Glist_q0, qgrid, ωgrid, syms, lattice_Å)
	elseif eps0folder isa AbstractString
		eps0folder = [eps0folder]
		eps0_q̂₀ = [CartesianCoordinates{Float64}(eps0_q̂₀)]
	end
	neps0folder = length(eps0folder)
	if eps0_ωindex == :default
		eps0_ωindex = fill!(ωindex, neps0folder)
	elseif any(x->length(x) ≠ nω, eps0_ωindex)
		error("Mismatched length of ωindex of eps and eps0!")
	end

	# matrices_q0, q0, Glist_q0
	q0 = Vector{CartesianCoordinates{Float64}}(undef, neps0folder + 1)
	q0[1] = CartesianCoordinates{Float64}(q̂₀)
	if isnothing(eps0_q̂₀)
		for ieps0folder in 1:neps0folder
			q0[ieps0folder+1] = q0[1]
		end
	end
	matrices_q0 = Matrix{Matrix{eltype(mat)}}(undef, neps0folder + 1, nω)
	matrices_q0[1, :] = matrices[q0_idx, :]
	Glist_q0 = Vector{Vector{ReducedCoordinates{Int}}}(undef, neps0folder + 1)
	Glist_q0[1] = Glist
	for ieps0folder in 1:neps0folder
		eps_head = NCDataset(joinpath(eps0folder[ieps0folder], "ndb." * suffix), "r"; format)
		Glist = eps_head["X_RL_vecs"][:, :]
		nG = size(Glist, 1)
		Glist = map(Base.OneTo(nG)) do iG
			round.(Int, rlattice \ CartesianCoordinates{Float64}(rlattice_parameter .* Glist[iG, :]))
		end
		Glist_q0[ieps0folder+1] = Glist
		eps = NCDataset(joinpath(epsfolder, "ndb." * suffix * "_fragment_$(q0_idx)"), "r"; format)
		matrices_iq = eps["X_Q_1"][:, :, :, eps0_ωindex[ieps0folder]]
		close(eps)
		for iω in Base.OneTo(nω)
			#note: yambo存储的是不包含单位矩阵部分，且完全对称化的介电函数。
			matrices_q0[ieps0folder+1, iω] = complex.(matrices_iq[1, :, :, iω], matrices_iq[2, :, :, iω])
		end
	end

	# use a same Glist for all q0.
	Glist_q0_unique = intersect(Glist_q0...)
	Glist_q0_unique = sort!(Glist_q0_unique; by = G->G ⋅ (rldot * G))
	for ieps0mat in 1:(neps0mat+1)
		G2idx = Dict(G => iG for (iG, G) in enumerate(Glist_q0[ieps0mat]))
		Gidx = map(G -> G2idx[G], Glist_q0_unique)
		for iω in Base.OneTo(nω)
			matrices_q0[ieps0mat, iω] = matrices_q0[ieps0mat, iω][Gidx, Gidx]
		end
	end

	return DielectricPWBasis(matrices_qirre, qirre, Glist_qirre, matrices_q0, q0, Glist_q0, qgrid, ωgrid, syms, lattice_Å)
end



