"""
	DirectLR{<:LongRangePart} <: DirectInterAction

Direct coulomb term between wannier basis only with long-range correction.
Fields:
- `norb::Int`: the number of orbitals;
- `LR::LongRangePart`: long-range part of direct term.
"""
struct DirectLR{L <: LongRangePart} <: DirectInterAction
	norb::Int
	LR::L
end
function (v::DirectLR)(k::ReducedCoordinates; isΓ::Bool = false, nearΓ::Bool = isΓ)
	A = zeros(ComplexF64, v.norb, v.norb)
	v.LR(Val(+), A, k; isΓ, nearΓ)
	return A
end
function (v::DirectLR)(k::ReducedCoordinates,
	nkorkgrid::Union{Integer, AbstractVector{<:Integer}, NTuple{3, <:Integer}}; isΓ::Bool = false, nearΓ::Bool = isΓ)
	A = zeros(ComplexF64, v.norb, v.norb)
	v.LR(Val(+), A, k, nkorkgrid; isΓ, nearΓ)
	return A
end
function (v::DirectLR)(A::AbstractMatrix, k::ReducedCoordinates; isΓ::Bool = false, nearΓ::Bool = isΓ)
	size(A) == (v.norb, v.norb) || error("Buffer size mismatch.")
	fill!(A, zero(ComplexF64))
	v.LR(Val(+), A, k; isΓ, nearΓ)
	return A
end
function (v::DirectLR)(A::AbstractMatrix, k::ReducedCoordinates,
	nkorkgrid::Union{Integer, AbstractVector{<:Integer}, NTuple{3, <:Integer}}; isΓ::Bool = false, nearΓ::Bool = isΓ)
	size(A) == (v.norb, v.norb) || error("Buffer size mismatch.")
	fill!(A, zero(ComplexF64))
	v.LR(Val(+), A, k, nkorkgrid; isΓ, nearΓ)
	return A
end
function (v::DirectLR)(::Val{:G0scale}, k::ReducedCoordinates)
	A = Matrix{ComplexF64}(undef, v.norb, v.norb)
	v.LR(Val(:G0scale), A, k)
	return A
end
function (v::DirectLR)(::Val{:G0scale}, A::AbstractMatrix, k::ReducedCoordinates)
	size(A) == (v.norb, v.norb) || error("Buffer size mismatch.")
	v.LR(Val(:G0scale), A, k)
	return A
end
(v::DirectLR)(::Val{:taylor_0}) = v.LR(Val(:taylor_0))
function Base.show(io::IO, U::DirectLR{L}) where {L}
	println(io, "Direct term only with long-range correction between electrionic wannier basis.")
end
