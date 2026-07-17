# A symmetry operation (SymOp) is a couple (W, w) of a matrix W and a translation w, in
# reduced coordinates, such that for each atom of type A at position a, W a + w is also an
# atom of type A.
# The matrix W is unitary in Cartesian coordinates, but not in reduced coordinates.
# This induces an operator in Reciprocal space
#   S = inv(W')

export SymOp, check_group, LattSymOp

struct ConjugateOperator end
Base.:(*)(::ConjugateOperator, c::Number) = conj(c)
Base.:(*)(c::Number, ::ConjugateOperator) = conj(c)
struct TimeReversalOperator end


struct SymOp
	W::Mat3{Int}
	w::Vec3{Float64}
	# in Reciprocal space
	S::Mat3{Int}
	# for time reversal, 1 represents no time-reversal and -1 represents time-reversal.
	t::Int
	#  旋量空间的转动矩阵
end
function SymOp(W::AbstractMatrix{<:Integer}, w::AbstractVector{<:Real}; t = 1)
	w = mod.(w, 1)
	# Accatually, it should be inv(W'), but in this programme we will apply it to the state, so we ommmit the inv.
	# After testing, W' and inv(W') both are ok.
	S = t * transpose(W)
	return SymOp(W, w, S, t)
end



# function Base.convert(::Type{SymOp}, S)
# 	return SymOp(S.W, S.w, S.S)
# end

Base.:(==)(op1::SymOp, op2::SymOp) = op1.W == op2.W && op1.w == op2.w
function Base.isapprox(op1::SymOp, op2::SymOp; symtol = 1e-5)
	return op1.W == op2.W && is_approx_integer(op1.w - op2.w; atol = symtol)
end

Base.one(::Type{SymOp}) = SymOp(Mat3{Int}(I), Vec3(zeros(Int, 3)))
Base.one(::SymOp) = one(SymOp)
Base.isone(op::SymOp) = isone(op.W) && iszero(op.w)


# group composition and inverse.
function Base.:(*)(op1::SymOp, op2::SymOp)
	W = op1.W * op2.W
	w = op1.W * op2.w + op1.w
	return SymOp(W, w)
end
function Base.:(*)(op::SymOp, vec::AbstractVector{<:Real})
	return op.W * vec + op.w
end
Base.inv(op::SymOp) = SymOp(inv(op.W), -op.W \ op.w)


function check_group(symops::AbstractVector{SymOp}; kwargs...)
	is_approx_in_symops(s1) = any(s -> isapprox(s, s1; kwargs...), symops)
	is_approx_in_symops(one(SymOp)) || error("check_group: no identity element")
	for s in symops
		if !is_approx_in_symops(inv(s))
			error("check_group: symop $s with inverse $(inv(s)) is not in the group")
		end
		for s′ in symops
			if !is_approx_in_symops(s * s′)
				error("check_group: product is not stable")
			end
		end
	end
	return symops
end


"""
	qmap_idx, qmap_self, qmap_sym, qmap_Gwind = _qre_qirre_map(qre, qirre, syms)

-`qre`::Vector{ReducedCoordinates{Rational{Int}}}
-`qirre`::Vector{ReducedCoordinates{Rational{Int}}}
-`syms`::Vector{SymOp}
Used to map irreducible qpoints to reducible qpoints with symmetry operations(use sym.S).
Only requires `qre` don't have repetitive elements.

If no qpoints in `qre` mapped by a `qirre`, then qmap_idx[iq] = 0.
"""
function _qre_qirre_map(qre, qirre, syms)

	nqre = length(qre)
	# qre[iq] = sym * qirre[qmap_idx[iq]] - Gwind[iq]
	qmap_idx = zeros(Int, nqre)
	qmap_self = Vector{Bool}(undef, nqre)
	qmap_sym = Vector{SymOp}(undef, nqre)
	qmap_Gwind = Vector{ReducedCoordinates{Int}}(undef, nqre)

	for iqirre in eachindex(qirre)
		q_from = qirre[iqirre]
		# if iszero(q_from)
		# 	q0idx = findfirst(q->all(isinteger, q), qre)
		# 	qmap_idx[q0idx] = iqirre
		# 	qmap_self[q0idx] = true
		# 	qmap_sym[q0idx] = SymOp([1 0 0; 0 1 0; 0 0 1], [0, 0, 0])
		# 	qmap_Gwind[q0idx] = ReducedCoordinates{Int}(-qre[q0idx])
		# 	continue
		# end
		for sym in syms
			q_from_sym = sym.S * q_from
			qidx = findfirst(qre) do q
				all(isinteger, q - q_from_sym)
			end
			# isnothing(qidx) && error("no qpoints in qre are mapped to ", q_from, "with ", sym)
			# 一些特殊的网格例如包含1//2偏移的网格会导致k网格在对称操作下不封闭。
			# 介电函数一般使用包含Γ点的k网格，但是保留该函数的通用性以及灵活性，故不封闭也不报错。
			isnothing(qidx) && continue
			if iszero(qmap_idx[qidx])
				qmap_idx[qidx] = iqirre
				if all(iszero, qre[qidx] - q_from)
					qmap_self[qidx] = true
					qmap_sym[qidx] = SymOp([1 0 0; 0 1 0; 0 0 1], [0, 0, 0])
					qmap_Gwind[qidx] = ReducedCoordinates{Int}(0, 0, 0)
				else
					qmap_self[qidx] = false
					qmap_sym[qidx] = sym
					qmap_Gwind[qidx] = q_from_sym - qre[qidx]
				end
			end
		end
	end

	return qmap_idx, qmap_self, qmap_sym, qmap_Gwind
end
"""
	qmap_idx, qmap_self, qmap_sym, qmap_Gwind = _qre_qirre_map!(qre, qirre, syms)

The same as `_qre_qirre_map!`.
Change `qre` if all(isinteger, qre - qirre).
"""
function _qre_qirre_map!(qre, qirre, syms)

	nqre = length(qre)
	# qre[iq] = sym * qirre[qmap_idx[iq]] - Gwind[iq]
	qmap_idx = zeros(Int, nqre)
	qmap_self = Vector{Bool}(undef, nqre)
	qmap_sym = Vector{SymOp}(undef, nqre)
	qmap_Gwind = Vector{ReducedCoordinates{Int}}(undef, nqre)

	for iqirre in eachindex(qirre)
		q_from = qirre[iqirre]
		for sym in syms
			q_from_sym = sym.S * q_from
			qidx = findfirst(qre) do q
				all(isinteger, q - q_from_sym)
			end
			# isnothing(qidx) && error("no qpoints in qre are mapped to ", q_from, "with ", sym)
			# 一些特殊的网格例如包含1//2偏移的网格会导致k网格在对称操作下不封闭。
			# 介电函数一般使用包含Γ点的k网格，但是保留该函数的通用性以及灵活性，故不封闭也不报错。
			isnothing(qidx) && continue
			if iszero(qmap_idx[qidx])
				qmap_idx[qidx] = iqirre
				if all(isinteger, qre[qidx] - q_from)
					qre[qidx] = q_from
					qmap_self[qidx] = true
					qmap_sym[qidx] = SymOp([1 0 0; 0 1 0; 0 0 1], [0, 0, 0])
					qmap_Gwind[qidx] = ReducedCoordinates{Int}(0, 0, 0)
				else
					qmap_self[qidx] = false
					qmap_sym[qidx] = sym
					qmap_Gwind[qidx] = q_from_sym - qre[qidx]
				end
			end
		end
	end

	return qmap_idx, qmap_self, qmap_sym, qmap_Gwind
end


# symop.S*k = k′, H(k′) = U(k)H(k)U(k)', |nk′⟩ = U(k)|nk⟩
struct LattSymOp{UT <: Number}
	k::ReducedCoordinates{Rational{Int}}
	k′::ReducedCoordinates{Rational{Int}}
	symop::SymOp
	U::Matrix{UT}
end

"""
Only work for one case presently: 
All orbitals are the same, and only one orbital(or two orbitals with opposite spin) on every site. 
In this case, we can ignore the factor of phase generated by rotating.
"""
function symop2lattsymop(symop, orblocation, k, k′; symtol = 1e-5)
	norb = length(orblocation)
	U = zeros(ComplexF64, norb, norb)
	for i in 1:norb
		aim = findfirst(x -> _is_approx_integer(symop * orblocation[i] - x; atol = symtol), orblocation)
		if isnothing(aim)
			error("Can't identify operated orbital.")
		end
		T = round.(symop * orblocation[i] - orblocation[aim])
		U[i, aim] = cispi(2 * k ⋅ T)
	end
	return LattSymOp(k, k′, symop, U)
end
