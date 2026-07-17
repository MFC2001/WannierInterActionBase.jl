export DielectricPWBasis, eps0sym
struct DielectricPWBasis{ET <: Number} <: AbstractDielectric
	# basic
	qgrid::RedKgrid # nq_re, q_re
	ωgrid::Vector{ComplexF64}
	nG::Vector{Int} # length.(Gidx)
	nG_max::Int # maximum(nG)
	Gidx::Vector{Vector{Int}} # [iq][nG[iq]].
	minGidx::Vector{Int} # nq_re.
	Glist::Vector{ReducedCoordinates{Int}} # G = Glist[Gidx]
	q0_idx::Int # findfirst(iszero, qgrid)
	lattice::Lattice{Float64}
	rlattice::ReciprocalLattice{Float64}
	ldot::Mat3{Float64}
	rldot::Mat3{Float64}
	syms::Vector{SymOp}
	# core data
	# matrices = ϵ⁻¹_{G₁,G₂}(q, ω) - δ_{G₁,G₂} = √V(q+G₁) χ_{G₁,G₂}(q, ω) √V(q+G₂)
	# only for q ≠ 0.
	matrices::Matrix{Matrix{ET}} # [nq_irre, nω]
	q_mat::Vector{ReducedCoordinates{Rational{Int}}} # nq_irre, q_irre
	qmap_idx::Vector{Int} # q_re -> q_irre; qmap_idx[q0_idx] = 0
	qmap_self::Vector{Bool} # q_re == q_irre
	qmap_sym::Vector{SymOp} # q_re[iq] = qmap_sym[iq].S * q_irre[qmap_idx[iq]] - qmap_Gwind[iq]
	qmap_Gwind::Vector{ReducedCoordinates{Int}}
	# for q = 0, 要求这些q0对应的介电矩阵G基需完全相同。
	matrices_q0::Matrix{Matrix{ET}} # [nq0, nω]
	q₀::Vector{CartesianCoordinates{Float64}}
	q̂₀::Vector{CartesianCoordinates{Float64}}
end
function Base.getproperty(eps::DielectricPWBasis, name::Symbol)
	if name === :nq
		return length(eps.qgrid)
	elseif name === :nω
		return length(eps.ωgrid)
	else
		return getfield(eps, name)
	end
end
function Base.show(io::IO, eps::DielectricPWBasis)
	print(io, "DielectricPWBasis with nq = $(eps.nq), nω = $(eps.nω)")
	print(io, ", nG_max = $(eps.nG_max)")
	return nothing
end
function Base.getindex(eps::DielectricPWBasis, iq::Int, ::Colon)
	return eps[iq, 1:eps.nω]
end
function Base.getindex(eps::DielectricPWBasis{ET}, iq::Int, iω::UnitRange{Int}) where {ET}
	iq == eps.q0_idx && return eps.matrices_q0[1, iω]
	eps.qmap_self[iq] && return eps.matrices[eps.qmap_idx[iq], iω]
	T = eps.qmap_sym[iq].w
	iszero(T) && return eps.matrices[eps.qmap_idx[iq], iω]
	# need to correct or optimize
	nG_q = eps.nG[iq]
	Glist_q = @view eps.Glist[eps.Gidx[iq]]
	mat_qirre = eps.matrices[eps.qmap_idx[iq], iω]
	mat_q = similar(mat_qirre)
	n_aimω = length(iω)
	for iω in Base.OneTo(n_aimω)
		mat_q[iω] = Matrix{ET}(undef, nG_q, nG_q)
	end
	for iG2 in 2:nG_q, iG1 in 1:(iG2-1)
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
end
function Base.getindex(eps::DielectricPWBasis{ET}, iq::Int, iω::Int) where {ET}
	iq == eps.q0_idx && return eps.matrices_q0[1, iω]
	eps.qmap_self[iq] && return eps.matrices[eps.qmap_idx[iq], iω]
	T = eps.qmap_sym[iq].w
	iszero(T) && return eps.matrices[eps.qmap_idx[iq], iω]
	#need to correct or optimize
	nG_q = eps.nG[iq]
	Glist_q = @view eps.Glist[eps.Gidx[iq]]
	mat_qirre = eps.matrices[eps.qmap_idx[iq], iω]
	mat_q = similar(mat_qirre)
	for iG2 in 2:nG_q, iG1 in 1:(iG2-1)
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
end
function Base.getindex(eps::DielectricPWBasis{ET}, idx::Int) where {ET}
	if 0 < idx < eps.nq * eps.nω
		iw = (idx - 1) ÷ eps.nq + 1
		iq = (idx - 1) % eps.nq + 1
		return eps[iq, iw]
	else
		error("Wrong idx for eps.")
	end
end
function (eps::DielectricPWBasis{ET})(::Val{:head}, iq::Int, iω::Int) where {ET}
	if iq == eps.q0_idx
		G0 = eps.G0idx[iq]
		return eps.matrices_q0[1, iω][G0, G0]
	else
		iq_irre = eps.qmap_idx[iq]
		G0 = eps.G0idx[iq]
		return eps.matrices[iq_irre, iω][G0, G0]
	end
end
function (eps::DielectricPWBasis{ET})(::Val{:head}, iq::Int, iω::UnitRange{Int}) where {ET}
	if iq == eps.q0_idx
		G0 = eps.G0idx[iq]
		return map(i->eps.matrices_q0[1, i][G0, G0], iω)
	else
		iq_irre = eps.qmap_idx[iq]
		G0 = eps.G0idx[iq]
		return map(i->eps.matrices[iq_irre, i][G0, G0], iω)
	end
end
function (eps::DielectricPWBasis{ET})(::Val{:head}, iq::Int, ::Colon) where {ET}
	return eps(Val(:head), iq, 1:eps.nω)
end
function eps0sym(eps::DielectricPWBasis{ET}) where {ET}
	#TODO
	# 见BerkeleyGW.
	error("To be continued!")
end
"""

-`matrices`::Matrix{Matrix{<:Number}}
-`qirre`::Vector{ReducedCoordinates{Rational{Int}}}
-`Glist_qirre`::Vector{Vector{ReducedCoordinates{Int}}}
-`matrices_q0`::Matrix{Matrix{<:Number}}
-`q₀`::Vector{CartesianCoordinates{Float64}}
-`Glist_q0`::Vector{ReducedCoordinates{Int}}
-`qgrid`::RedKgrid
-`ωgrid`::Union{Number, Vector{<:Number}}
-`syms`::Vector{SymOp}

-`lattice`::Lattice{<:Real}
"""
function DielectricPWBasis(matrices, qirre, Glist_qirre, matrices_q0, q₀, Glist_q0, qgrid, ωgrid, syms, lattice)
	if ωgrid isa Number
		ωgrid = [ωgrid]
	end
	nω = length(ωgrid)

	rlattice = reciprocal(lattice)
	ldot = transpose(parent(lattice)) * parent(lattice)
	rldot = transpose(parent(rlattice)) * parent(rlattice)

	(qmap_idx, qmap_self, qmap_sym, qmap_Gwind) = _qre_qirre_map(qgrid.kdirect, qirre, syms)

	# qgrid should be Γ-centered grid.
	q0_idx = findfirst(iszero, qgrid)
	isnothing(q0_idx) && error("No Γ in qgrid!")
	nomap_idx = findall(iszero, qmap_idx)
	if !(length(nomap_idx) == 1 && nomap_idx[1] == q0_idx)
		error("Wrong qpoints or syms!")
	end

	(Glist, Gidx) = _Gidx_irre2grid_withq0(Glist_qirre, Glist_q0, qmap_idx, qmap_self, qmap_sym, qmap_Gwind, rldot)
	nG = length.(Gidx)
	nG_max = maximum(nG)
	minGidx = map(Gidx) do iGs
		# Glist 已经按照模长排序过。
		(_, minGidx_iq) = findmin(iGs)
		return minGidx_iq
	end

	nq₀ = length(q₀)
	q₀_new = Vector{CartesianCoordinates{Float64}}(undef, nq₀)
	q̂₀ = Vector{CartesianCoordinates{Float64}}(undef, nq₀)
	for iq₀ in 1:nq₀
		normq₀ = sqrt(q₀[iq₀] ⋅ q₀[iq₀])
		if isapprox(normq₀, 1; atol = 1e-2)
			# q₀[iq₀] is a unit vector.
			q₀_new[iq₀] = CartesianCoordinates{Float64}(0, 0, 0)
			q̂₀[iq₀] = CartesianCoordinates{Float64}(q₀[iq₀] ./ normq₀)
		elseif isapprox(normq₀, 0; atol = 1e-8)
			# q₀[iq₀] is a zero vector.
			# q₀_new[iq₀] = CartesianCoordinates{Float64}(q₀[iq₀])
			# q̂₀[iq₀] = CartesianCoordinates{Float64}(0, 0, 0)
			error("Don't accept a zero vector as q₀!")
		else
			# q₀[iq₀] is a small vector.
			q₀_new[iq₀] = CartesianCoordinates{Float64}(q₀[iq₀])
			q̂₀[iq₀] = CartesianCoordinates{Float64}(q₀[iq₀] ./ normq₀)
		end
	end

	return DielectricPWBasis(qgrid, ωgrid,
		nG, nG_max, Gidx, minGidx, Glist,
		q0_idx, lattice, rlattice, ldot, rldot, syms,
		matrices, qirre,
		qmap_idx, qmap_self, qmap_sym, qmap_Gwind,
		matrices_q0, q₀_new, q̂₀,
	)
end
function _Gidx_irre2grid_withq0(Glist_qirre, Glist_q0, qmap_idx, qmap_self, qmap_sym, qmap_Gwind, rldot)
	nqre = length(qmap_idx)
	#Glist_q
	Glist_qgrid = Vector{Vector{ReducedCoordinates{Int}}}(undef, nqre)
	for (iq_re, iq_irre) in enumerate(qmap_idx)
		if iszero(iq_irre)
			Glist_qgrid[iq_re] = Glist_q0
		else
			if qmap_self[iq_re]
				Glist_qgrid[iq_re] = Glist_qirre[iq_irre]
			else
				Glist_qgrid[iq_re] = map(Glist_qirre[iq_irre]) do G
					qmap_sym[iq_re].S * G + qmap_Gwind[iq_re]
				end
			end
		end
	end

	# merge Glist
	Glist_qgrid_unique = Set(Iterators.flatten(Glist_qgrid))
	Glist_qgrid_unique = collect(Glist_qgrid_unique)
	sort!(Glist_qgrid_unique; by = G -> G ⋅ (rldot * G))

	Glist_qgrid_unique_Dict = Dict(G => i for (i, G) in enumerate(Glist_qgrid_unique))

	Gidx_qgrid = Vector{Vector{Int}}(undef, nqre)
	for iq in Base.OneTo(nqre)
		Gidx_qgrid[iq] = map(G -> Glist_qgrid_unique_Dict[G], Glist_qgrid[iq])
	end

	return Glist_qgrid_unique, Gidx_qgrid
end
