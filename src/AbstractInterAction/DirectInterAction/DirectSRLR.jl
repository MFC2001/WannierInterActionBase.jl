"""
	DirectSRLR{<: AbstractReciprocalHoppings, <: LongRangePart} <: DirectInterAction

Direct coulomb term between wannier basis with long-range correction.
Fields:
- `norb::Int`: the number of orbitals;
- `SR::AbstractReciprocalHoppings`: short-range part of direct term;
- `LR::LongRangePart`: long-range part of direct term.

We have achieved long-range correction by using Gauss potential.
If you want to use other potential, you can define a new subtype of `LongRangePart`. 
Make sure your subtype has methods:

	(v::yoursubtype)(::Val{+}, paras...)
	(v::yoursubtype)(::Val{+}, paras...)

For an instance `u` of `DirectSRLR`, you can run:

```julia
julia> u(k)
julia> u(k, nk)
```

These two methods will create a new matrix. And `u(k, nk)` will judge whether head term or normal value by k.

```julia
julia> u(A, k)
julia> u(A, k, nk)
```

These two methods will change `A` instead of creating a new matrix.
"""
struct DirectSRLR{S <: AbstractReciprocalHoppings, L <: LongRangePart} <: DirectInterAction
	norb::Int
	SR::S
	LR::L
end
function (v::DirectSRLR)(k::ReducedCoordinates; isΓ::Bool = false, nearΓ::Bool = isΓ)
	A = Matrix{ComplexF64}(undef, v.norb, v.norb)
	return v(A, k; isΓ, nearΓ)
end
function (v::DirectSRLR)(k::ReducedCoordinates,
	nkorkgrid::Union{Integer, AbstractVector{<:Integer}, NTuple{3, <:Integer}}; isΓ::Bool = false, nearΓ::Bool = isΓ)
	A = Matrix{ComplexF64}(undef, v.norb, v.norb)
	return v(A, k, nkorkgrid; isΓ, nearΓ)
end
function (v::DirectSRLR)(A::AbstractMatrix, k::ReducedCoordinates; isΓ::Bool = false, nearΓ::Bool = isΓ)
	isΓ ? v.SR(A, ReducedCoordinates(0, 0, 0)) : v.SR(A, k)
	v.LR(Val(+), A, k; isΓ, nearΓ)
	return A
end
function (v::DirectSRLR)(A::AbstractMatrix, k::ReducedCoordinates,
	nkorkgrid::Union{Integer, AbstractVector{<:Integer}, NTuple{3, <:Integer}}; isΓ::Bool = false, nearΓ::Bool = isΓ)
	isΓ ? v.SR(A, ReducedCoordinates(0, 0, 0)) : v.SR(A, k)
	v.LR(Val(+), A, k, nkorkgrid; isΓ, nearΓ)
	return A
end
function (v::DirectSRLR)(::Val{:G0scale}, k::ReducedCoordinates)
	A = Matrix{ComplexF64}(undef, v.norb, v.norb)
	v.LR(Val(:G0scale), A, k)
	return A
end
function (v::DirectSRLR)(::Val{:G0scale}, A::AbstractMatrix, k::ReducedCoordinates)
	size(A) == (v.norb, v.norb) || error("Buffer size mismatch.")
	v.LR(Val(:G0scale), A, k)
	return A
end
(v::DirectSRLR)(::Val{:taylor_0}) = v.LR(Val(:taylor_0))
function Base.show(io::IO, U::DirectSRLR{S, L}) where {S, L}
	println(io, "Direct term with long-range correction between electrionic wannier basis.")
end
function _DirectSRLR_SR_3D(U, φ, lattice, orblocat)

	SR = deepcopy(U)
	hops_SR = hops(SR)

	for I in CartesianIndices(hops_SR)
		(i, j) = Tuple(I)
		dorb = orblocat[j] - orblocat[i]
		φᵢⱼ = φ[i, j]
		for (ii, hop) in enumerate(hops_SR[I])
			value_SR = value(hop) - φᵢⱼ(lattice * (path(hop) + dorb))
			hops_SR[I][ii] = similar(hop, value_SR)
		end
	end

	return SR
end
function _DirectSRLR_SR_2D(U, φ, lattice, orblocat, xyindex)

	lattice2D = SMatrix{2, 2, Float64, 4}(lattice[xyindex, xyindex])

	SR = deepcopy(U)
	hops_SR = hops(SR)

	for I in CartesianIndices(hops_SR)
		(i, j) = Tuple(I)
		dorb = orblocat[j] - orblocat[i]
		φᵢⱼ = φ[i, j]
		for (ii, hop) in enumerate(hops_SR[I])
			rfrac = path(hop) + dorb
			rfrac = rfrac[xyindex]
			value_SR = value(hop) - φᵢⱼ(lattice2D * rfrac)
			hops_SR[I][ii] = similar(hop, value_SR)
		end
	end

	return SR
end
function DirectSRLR(::Val{:gauss}, U::AbstractReciprocalHoppings,
	lattice::Lattice, orblocat::Vector{<:ReducedCoordinates}; norb::Integer = length(orblocat),
	α::Real = 1, αmat::AbstractMatrix{<:Real} = fill(α, (norb, norb)), δ::Real = 1e-8, ϵ::Real = 1)

	αmat = _LongRangePartGauss_symmetrize_αmat(αmat)
	φR = [RealCoulomb(:gauss; ϵ, α = αmat[i, j]) for i in 1:norb, j in 1:norb]
	V_SR = _DirectSRLR_SR_3D(U, φR, lattice, orblocat)
	V_LR = LongRangePart(Val(:gauss), lattice, orblocat; norb, αmat, δ, ϵ)

	return DirectSRLR(norb, V_SR, V_LR)
end
function DirectSRLR(::Val{:keldysh_gauss}, U::AbstractReciprocalHoppings,
	lattice::Lattice, orblocat::Vector{<:ReducedCoordinates}; norb::Integer = length(orblocat),
	period = Bool[1, 1, 0],
	α::Real = 1, αmat::AbstractMatrix{<:Real} = fill(α, (norb, norb)),
	αz::Real = 1, αzmat::AbstractMatrix{<:Real} = fill(αz, (norb, norb)),
	r₀::Real = 0, ϵ₁::Real = 1, ϵ₂::Real = 1, r₀atol::Real = 0, δ::Real = 1e-8)

	period = Vec3{Bool}(period)
	if !period[1]
		xyindex = SVector{2, Int}(2, 3)
		zindex = 1
	elseif !period[2]
		xyindex = SVector{2, Int}(3, 1)
		zindex = 2
	else
		xyindex = SVector{2, Int}(1, 2)
		zindex = 3
	end

	αmat, αzmat = _LongRangePartKeldyshGauss_symmetrize_αmat(αmat, αzmat)
	φR = [RealCoulomb(:keldysh_gauss; r₀, ϵ₁, ϵ₂, α = αmat[i, j], αz = αzmat[i, j], atol = r₀atol) for i in 1:norb, j in 1:norb]
	V_SR = _DirectSRLR_SR_2D(U, φR, lattice, orblocat, xyindex)
	V_LR = LongRangePart(Val(:keldysh_gauss), lattice, orblocat; norb, αmat, αzmat, r₀, ϵ₁, ϵ₂, r₀atol, δ)

	return DirectSRLR(norb, V_SR, V_LR)
end
