export WanIntPairType, WanIntPairVecType, WanIntPairIdxType,
	atom2orbpair, typeofpair, completepair, uplopairsym, uplopairidx

# Used for function definition.
const WanIntPairVecType = Tuple{<:Integer, <:Integer, <:AbstractVector{<:Integer}}
const WanIntPairIdxType = Tuple{<:Integer, <:Integer, <:Integer}
const WanIntPairType = Union{WanIntPairVecType, WanIntPairIdxType}

"""
	atom2orbpair(cell::Cell, rcut::Real, atom_orb_index; shell_tol::Real = 1e-6)

Generate `atompair` firstly, then use atom2orbpair(atompair, atom_orb_index)

	atom2orbpair(atompair, atom_orb_index)

# Examples
```julia
julia> atompair = [(1, 1, 0, 0, 0)];
julia> atom_orb_index = [[1, 2]];
julia> atom2orbpair(atompair, atom_orb_index)
4-element Vector{NTuple{5, Int64}}:
 (1, 1, 0, 0, 0)
 (2, 1, 0, 0, 0)
 (1, 2, 0, 0, 0)
 (2, 2, 0, 0, 0)
```
"""
function atom2orbpair(cell::Cell, rcut::Real, atom_orb_index; shell_tol::Real = 1e-6)
	# Generate atom pairs firstly.
	neighborlist, neighbordist = findneighbors(cell, rcut; shell_tol);
	atom_pair = Iterators.flatten(Iterators.flatten(neighborlist))
	# Obtain orbital pairs from atomic pairs.
	return atom2orbpair(atom_pair, atom_orb_index)
end
function atom2orbpair(atompair, atom_orb_index)
	orbpair = Vector{NTuple{5, Int}}(undef, 0)
	for p in atompair
		pp = vec([(i, j, p[3], p[4], p[5]) for i in atom_orb_index[p[1]], j in atom_orb_index[p[2]]])
		append!(orbpair, pp)
	end
	return orbpair
end
"""
	typeofpair(pair)
	typeofpair(pair, pR)
	typeofpair(pairs)
	typeofpair(pairs, pR)

Obtain the type of pair: :0, :p, :n.
"""
typeofpair(pair::AbstractVector{<:WanIntPairVecType}) = map(p->typeofpair(p), pair)
typeofpair(pair::AbstractVector{<:WanIntPairIdxType}, pR) = map(p->typeofpair(p, pR), pair)
function typeofpair(pair::WanIntPairVecType)::Symbol
	_iszero_pair(pair) && return Symbol(0)
	_ispositive_pair(pair) && return :p
	_isnegative_pair(pair) && return :n
	error("Can't define the type of $pair")
end
function typeofpair(pair::WanIntPairIdxType, pR::AbstractVector{<:AbstractVector{<:Integer}})::Symbol
	_iszero_pair(pair, pR) && return Symbol(0)
	_ispositive_pair(pair, pR) && return :p
	_isnegative_pair(pair, pR) && return :n
	error("Can't define the type of $pair")
end
"""
	completepair(pairs)
	completepair(pairs, pR)

Obtain the type of pair: :0, :p, :n.
"""
function completepair(pair::AbstractVector{<:WanIntPairVecType})
	original_pair = [(Int(p[1]), Int(p[2]), ReducedCoordinates{Int}(p[3])) for p in pair]
	conj_pair = [(p[2], p[1], -p[3]) for p in original_pair]
	by(p) = (norm2(p[3]), p[1], p[2])
	pair = sort(unique([original_pair; conj_pair]); by)
	return pair
end
function completepair(pair::AbstractVector{<:WanIntPairIdxType}, pR::AbstractVector{<:AbstractVector{<:Integer}})
	# obtain pR and conjpR, and original_pR_idx
	original_pR = ReducedCoordinates{Int}.(collect(pR))
	conj_pR = [-RR for RR in original_pR]
	pR = sort(unique([original_pR, conj_pR]); by = norm2)
	pR2idx = Dict(x=>i for (i, x) in enumerate(pR))
	original_pR_idx = [pR2idx[RR] for RR in original_pR]
	conjpR = [pR2idx[-RR] for RR in pR]
	# obtain pair
	original_pair = [(Int(p[1]), Int(p[2]), original_pR_idx[p[3]]) for p in pair]
	conj_pair = [(p[2], p[1], conjpR[p[3]]) for p in original_pair]
	by(p) = (p[3], p[1], p[2])
	pair = sort(unique([original_pair; conj_pair]); by)
	return pair, pR
end
"""
	uplopairsym(uplo)
"""
uplopairsym(uplo::Symbol) = uplopairsym(Val(uplo))
uplopairsym(::Val{:U}) = Set([(Symbol(0), :p), (:p, Symbol(0)), (:p, :p), (:p, :n)])
uplopairsym(::Val{:L}) = Set([(Symbol(0), :n), (:n, Symbol(0)), (:n, :n), (:n, :p)])
uplopairsym(::Val{:all}) = union!(uplopairsym(Val(:U)), uplopairsym(Val(:L)))
"""
	uplopairidx(pair, uplo = :U)
	uplopairidx(pair, pR, uplo = :U)
	uplopairidx(pairtype, uplo = :U)
"""
function uplopairidx(pair::AbstractVector{<:WanIntPairVecType}, uplo::Symbol = :U)
	pairtype = map(p->typeofpair(p), pair)
	return uplopairidx(pairtype, uplo)
end
function uplopairidx(pair::AbstractVector{<:WanIntPairIdxType}, pR::AbstractVector{<:AbstractVector{<:Integer}}, uplo::Symbol = :U)
	pairtype = map(p->typeofpair(p, pR), pair)
	return uplopairidx(pairtype, uplo)
end
function uplopairidx(pairtype::AbstractVector{Symbol}, uplo::Symbol = :U)
	uplopair = uplopairsym(Val(uplo))
	np = length(pairtype)
	pair_idx = vec([(ip1, ip2) for ip1 in Base.OneTo(np), ip2 in Base.OneTo(np)])
	filter!(pair_idx) do (ip1, ip2)
		(pairtype[ip1], pairtype[ip2]) ∈ uplopair
	end
	return pair_idx
end

"""
	_iszero_pair(pair)
	_iszero_pair(pair, pR)
"""
@inline _iszero_pair(pair::WanIntPairVecType) = pair[1] == pair[2] && iszero(pair[3])
@inline _iszero_pair(pair::WanIntPairIdxType, pR::AbstractVector{<:AbstractVector{<:Integer}}) = pair[1] == pair[2] && iszero(pR[pair[3]])
"""
	_ispositive_pair(pair)
	_ispositive_pair(pair, pR)
"""
@inline function _ispositive_pair(pair::WanIntPairVecType)
	pair[1] < pair[2] && return true
	pair[1] == pair[2] && return _ispositive_vector(pair[3])
	return false
end
@inline function _ispositive_pair(pair::WanIntPairIdxType, pR::AbstractVector{<:AbstractVector{<:Integer}})
	pair[1] < pair[2] && return true
	pair[1] == pair[2] && return _ispositive_vector(pR[pair[3]])
	return false
end
"""
	_isnegative_pair(pair)
	_isnegative_pair(pair, pR)
"""
@inline function _isnegative_pair(pair::WanIntPairVecType)
	pair[1] > pair[2] && return true
	pair[1] == pair[2] && return _isnegative_vector(pair[3])
	return false
end
@inline function _isnegative_pair(pair::WanIntPairIdxType, pR::AbstractVector{<:AbstractVector{<:Integer}})
	pair[1] > pair[2] && return true
	pair[1] == pair[2] && return _isnegative_vector(pR[pair[3]])
	return false
end
"""
	_isvalid_pair(pair, nw)
	_isvalid_pair(pair, nw, npR)
"""
@inline _isvalid_pair(pair::WanIntPairType, nw::Integer) = 1 ≤ pair[1] ≤ nw && 1 ≤ pair[2] ≤ nw
@inline _isvalid_pair(pair::WanIntPairIdxType, nw::Integer, npR::Integer) = 1 ≤ pair[1] ≤ nw && 1 ≤ pair[2] ≤ nw && 1 ≤ pair[3] ≤ npR
