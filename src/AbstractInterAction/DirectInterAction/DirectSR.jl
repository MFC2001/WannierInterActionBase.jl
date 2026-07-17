"""
	DirectSR{<: AbstractReciprocalHoppings} <: DirectInterAction

Direct coulomb term between wannier basis without long-range correction.
Fields:
- `norb::Int`: the number of orbitals;
- `SR::AbstractReciprocalHoppings`: short-range part of direct term.
"""
struct DirectSR{S <: AbstractReciprocalHoppings} <: DirectInterAction
	norb::Int
	SR::S
end
function (v::DirectSR)(k::ReducedCoordinates; isΓ::Bool = false, nearΓ::Bool = isΓ)
	A = Matrix{ComplexF64}(undef, v.norb, v.norb)
	return v(A, k; isΓ, nearΓ)
end
function (v::DirectSR)(k::ReducedCoordinates,
	nkorkgrid::Union{Integer, AbstractVector{<:Integer}, NTuple{3, <:Integer}}; isΓ::Bool = false, nearΓ::Bool = isΓ)
	A = Matrix{ComplexF64}(undef, v.norb, v.norb)
	return v(A, k, nkorkgrid; isΓ, nearΓ)
end
(v::DirectSR)(A::AbstractMatrix, k::ReducedCoordinates; isΓ::Bool = false, nearΓ::Bool = isΓ) =
	isΓ ? v.SR(A, ReducedCoordinates(0, 0, 0)) : v.SR(A, k)
(v::DirectSR)(A::AbstractMatrix, k::ReducedCoordinates, nkorkgrid::Union{Integer, AbstractVector{<:Integer}, NTuple{3, <:Integer}}; isΓ::Bool = false, nearΓ::Bool = isΓ) =
	isΓ ? v.SR(A, ReducedCoordinates(0, 0, 0)) : v.SR(A, k)
(v::DirectSR)(::Val{:G0scale}, k::ReducedCoordinates) = zeros(ComplexF64, v.norb, v.norb)
function (v::DirectSR)(::Val{:G0scale}, A::AbstractMatrix, k::ReducedCoordinates)
	size(A) == (v.norb, v.norb) || error("Buffer size mismatch.")
	fill!(A, zero(ComplexF64))
	return A
end
(v::DirectSR)(::Val{:taylor_0}) = zeros(Float64, v.norb, v.norb)
function Base.show(io::IO, U::DirectSR{S}) where {S}
	println(io, "Direct term without long-range correction between electrionic wannier basis.")
end
