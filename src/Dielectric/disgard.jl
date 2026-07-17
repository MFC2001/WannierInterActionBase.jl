struct DielectricMatrixNone <: DielectricMatrix
	qgrid::RedKgrid
	q0::ReducedCoordinates{Float64}
	q0_idx::Int
	have_negq0::Bool
end
function DielectricMatrixNone(; kgrid_size, q0 = ReducedCoordinates(0.0, 0.0, 1e-6))
	qgrid = RedKgrid(MonkhorstPack(kgrid_size))
	q0_idx = findfirst(iszero, qgrid)
	have_negq0 = false
	return DielectricMatrixNone(qgrid, q0, q0_idx, have_negq0)
end

struct DielectricMatrix_scalar{ET <: Number} <: DielectricMatrix
	#public
	nω::Int
	nq::Int
	qgrid::RedKgrid # nq_re, q_re
	ωgrid::Vector{ComplexF64}
	ϵ00::Matrix{ET} # [nq, nω]
	nG::Vector{Int} # length.(Gidx)
	nG_max::Int # maximum(nG)
	Gidx::Vector{Vector{Int}} # [iq][nG[iq]].
	Glist::Vector{ReducedCoordinates{Int}} # G = Glist[Gidx]
	q0::ReducedCoordinates{Float64}
	q0_idx::Int
	have_negq0::Bool
	nG_negq0::Int # length.(Gidx_negq0)
	Gidx_negq0::Vector{Int} # [nG_negq0]
	qG_norm::Vector{Vector{Float64}} # [iq][iG]
	qG_norm_negq0::Vector{Float64} # [iG_negq0]
	lattice::Lattice{Float64}
	rlattice::ReciprocalLattice{Float64}
	ldot::Mat3{Float64}
	rldot::Mat3{Float64}
	Ω::Float64
	rΩ::Float64
	#private
	#note: 这里存储的介电矩阵为epsilon - 1，即已经去掉了单位矩阵部分，仅包含了极化函数部分。
	#note: 这里存储的是完全对称化的介电矩阵，包括截断因子也开方对称分配。
	mat_inv::Matrix{Matrix{ET}} # [nq_irre, nω][nG[iq_irre], nG[iq_irre]]
	q_mat::Vector{ReducedCoordinates{Rational{Int}}} # nq_irre
	qmap_idx::Vector{Int} # q_re -> q_irre
	qmap_self::Vector{Bool} # q_re == q_irre
	qmap_sym::Vector{SymOp} # q_re[iq] = qmap_sym[iq] * q_irre[qmap_idx[iq]] - qmap_Gwind[iq]
	qmap_Gwind::Vector{ReducedCoordinates{Int}}
	mat_inv_negq0::Vector{Matrix{ET}} # [nω][nG_negq0, nG_negq0]
end
function Base.show(io::IO, eps::DielectricMatrix_scalar)
	println(io, "DielectricMatrix_scalar with nq = $(eps.nq), nω = $(eps.nω)")
	println(io, "nG_max = $(eps.nG_max)")
	eps.have_negq0 ? println(io, "with negq0") : println(io, "without negq0")
	return nothing
end
function Base.getindex(eps::DielectricMatrix_scalar{ET}, q::Int, ::Colon) where {ET}
	if q > 0
		return eps[q, 1:eps.nω]
	elseif q == -1
		if eps.have_negq0
			return eps.mat_inv_negq0
		else
			error("negq0 is not available")
		end
	else
		error("invalid q index")
	end
end
function Base.getindex(eps::DielectricMatrix_scalar{ET}, q::Int, ω::UnitRange{Int}) where {ET}

	if q > 0

		eps.qmap_self[q] && return eps.mat_inv[eps.qmap_idx[q], ω]

		#need to correct or optimize
		mat_qirre = eps.mat_inv[eps.qmap_idx[q], ω]
		T = eps.qmap_sym[q].w
		Glist_q = @view eps.Glist[eps.Gidx[q]]
		nG_q = eps.nG[q]
		mat_q = similar(mat_qirre)

		n_aimω = length(ω)
		for iω in Base.OneTo(n_aimω)
			mat_q[iω] = Matrix{ET}(undef, nG_q, nG_q)
		end

		pairs = Vector{Tuple{Int, Int}}(undef, Int((nG_q * (nG_q - 1)) // 2))
		n = 0
		for iG2 in 2:nG_q, iG1 in 1:(iG2-1)
			n += 1
			pairs[n] = (iG1, iG2)
		end

		Threads.@threads for (iG1, iG2) in pairs
			G1 = Glist_q[iG1]
			G2 = Glist_q[iG2]
			phase = cispi(2 * (G1 - G2) ⋅ T)
			for iω in Base.OneTo(n_aimω)
				mat_q[iω][iG1, iG2] = mat_qirre[iω][iG1, iG2] * phase
				mat_q[iω][iG2, iG1] = mat_qirre[iω][iG2, iG1] * conj(phase)
			end
		end
		for iω in Base.OneTo(n_aimω)
			for iG in 1:nG_q
				mat_q[iω][iG, iG] = mat_qirre[iω][iG, iG]
			end
		end

		return mat_q

	elseif q == -1
		if eps.have_negq0
			return eps.mat_inv_negq0[ω]
		else
			error("negq0 is not available")
		end
	else
		error("invalid q index")
	end

end
function Base.getindex(eps::DielectricMatrix_scalar, q::Int, ω::Int)

	if q > 0

		eps.qmap_self[q] && return eps.mat_inv[eps.qmap_idx[q], ω]

		#need to correct or optimize
		mat_qirre = eps.mat_inv[eps.qmap_idx[q], ω]
		T = eps.qmap_sym[q].w
		Glist_q = @view eps.Glist[eps.Gidx[q]]
		nG_q = eps.nG[q]
		mat_q = similar(mat_qirre)

		pairs = Vector{Tuple{Int, Int}}(undef, Int((nG_q * (nG_q - 1)) // 2))
		n = 0
		for iG2 in 2:nG_q, iG1 in 1:(iG2-1)
			n += 1
			pairs[n] = (iG1, iG2)
		end

		Threads.@threads for (iG1, iG2) in pairs
			G1 = Glist_q[iG1]
			G2 = Glist_q[iG2]
			phase = cispi(2 * (G1 - G2) ⋅ T)
			mat_q[iG1, iG2] = mat_qirre[iG1, iG2] * phase
			mat_q[iG2, iG1] = mat_qirre[iG2, iG1] * conj(phase)
		end
		for iG in 1:nG_q
			mat_q[iG, iG] = mat_qirre[iG, iG]
		end

		return mat_q

	elseif q == -1
		if eps.have_negq0
			return eps.mat_inv_negq0[ω]
		else
			error("negq0 is not available")
		end
	end

end
function Base.getindex(eps::DielectricMatrix_scalar, idx::Int)
	if idx > 0
		iw = (idx - 1) ÷ eps.nq + 1
		iq = (idx - 1) % eps.nq + 1
		return eps[iq, iw]
	else
		return eps.mat_inv_negq0[-idx]
	end
end
function eps0sym(eps::DielectricMatrix_scalar)
	#TODO
	error("To be continued!")
end
function DielectricMatrix_scalar(mat_qω_inv, qirre, qgrid, syms, Glist_qirre, ωgrid, lattice, q0, qG_norm; #use q0 to calculate qG_norm.
	nω = length(ωgrid),
	have_negq0 = false,
	Glist_negq0 = Vector{ReducedCoordinates{Int}}(undef, 0),
	qG_norm_negq0 = Vector{Float64}(undef, 0),
	mat_inv_negq0 = Vector{Matrix{ComplexF64}}(undef, nω),
)

	if ωgrid isa Number
		ωgrid = [ωgrid]
	end

	rlattice = reciprocal(lattice)
	ldot = transpose(parent(lattice)) * parent(lattice)
	rldot = transpose(parent(rlattice)) * parent(rlattice)
	Ω = det(parent(lattice))
	rΩ = det(parent(rlattice))

	nqirre = length(qirre)
	ϵ00_irre = Matrix{eltype(eltype(mat_qω_inv))}(undef, nqirre, nω)
	for iq in 1:nqirre
		G0idx_iq = findfirst(iszero, Glist_qirre[iq])
		if isnothing(G0idx_iq)
			error("Not find G=0 in Glist_qirre for iq = $iq when read epsmat.")
		end
		for iω in 1:nω
			ϵ00_irre[iq, iω] = mat_qω_inv[iq, iω][G0idx_iq, G0idx_iq]
		end
	end

	# qirre 中必须存在一个零矢量。
	q0_idx = findfirst(iszero, qgrid)
	(qgrid, qmap_idx, qmap_self, qmap_sym, qmap_Gwind) = _qirre_qgrid_map(qirre, qgrid, syms)
	qG_norm = qG_norm[qmap_idx]
	if have_negq0
		(Glist, Gidx_qgrid, Gidx_negq0) =
			_Gidx_irre2grid_withq0(Glist_qirre, qmap_idx, qmap_self, qmap_sym, qmap_Gwind, rldot; Glist_negq0)
		nG_qgrid = length.(Gidx_qgrid)
		nG_negq0 = length(Gidx_negq0)
		nG_qgrid_max = maximum([nG_qgrid; nG_negq0])
	else
		Gidx_negq0 = Vector{Int}(undef, 0)
		nG_negq0 = 0
		(Glist, Gidx_qgrid) = _Gidx_irre2grid_withq0(Glist_qirre, qmap_idx, qmap_self, qmap_sym, qmap_Gwind, rldot)
		nG_qgrid = length.(Gidx_qgrid)
		nG_qgrid_max = maximum(nG_qgrid)
	end

	nq = length(qgrid)
	ϵ00 = Matrix{eltype(ϵ00_irre)}(undef, nq, nω)
	for iq in Base.OneTo(nq)
		ϵ00[iq, :] .= ϵ00_irre[qmap_idx[iq], :]
	end

	return DielectricMatrix_scalar(
		nω, nq, qgrid, ωgrid, ϵ00,
		nG_qgrid, nG_qgrid_max, Gidx_qgrid, Glist,
		q0, q0_idx,
		have_negq0, nG_negq0, Gidx_negq0,
		qG_norm, qG_norm_negq0,
		lattice, rlattice,
		ldot, rldot,
		Ω, rΩ,
		mat_qω_inv, qirre,
		qmap_idx, qmap_self, qmap_sym, qmap_Gwind,
		mat_inv_negq0,
	)
end
