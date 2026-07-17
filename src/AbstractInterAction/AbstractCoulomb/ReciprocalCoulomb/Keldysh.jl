abstract type ReciprocalCoulombKeldysh <: ReciprocalCoulomb end
struct ReciprocalCoulombKeldyshn0 <: ReciprocalCoulombKeldysh
	r₀::Float64
	ϵ₁::Float64
	ϵ₂::Float64
	S::Float64
	coff::Float64
end
function (v::ReciprocalCoulombKeldyshn0)(k::AbstractVector{<:Real}; atol::Real = 1e-10)::Float64
	k = sqrt(k ⋅ k)
	return k ≤ atol ? 0.0 : v.coff / (1 + v.r₀ * k) / k
end
function (v::ReciprocalCoulombKeldyshn0)(k::Real; atol::Real = 1e-10)::Float64
	signbit(k) && error("k should be positive")
	return k ≤ atol ? 0.0 : v.coff / (1 + v.r₀ * k) / k
end
function (v::ReciprocalCoulombKeldyshn0)(::Val{:k²}, k²; atol::Real = 1e-10)::Float64
	signbit(k²) && error("k² should be positive")
	return k² ≤ atol ? 0.0 : begin
		k = sqrt(k²)
		v.coff / (1 + v.r₀ * k) / k
	end
end
function (v::ReciprocalCoulombKeldyshn0)(::Val{:kscale}, k::AbstractVector{<:Real})::Float64
	k = sqrt(k ⋅ k)
	return v.coff / (1 + v.r₀ * k)
end
function (v::ReciprocalCoulombKeldyshn0)(::Val{:kscale}, k::Real)::Float64
	signbit(k) && error("k should be positive")
	return v.coff / (1 + v.r₀ * k)
end
function (v::ReciprocalCoulombKeldyshn0)(::Val{:head}, nk::Integer)::Float64
	signbit(nk) && error("nk should be positive")
	NS = nk * v.S
	qsz = sqrt(4π / NS)
	return v.coff * NS * log(1 + v.r₀ * qsz) / (2 * π * v.r₀)
end
function (v::ReciprocalCoulombKeldyshn0)(::Val{:head}, rlattice::ReciprocalLattice, kgrid::Union{AbstractVector, Tuple},
	zidx::Integer = findfirst(nk->nk == 1, kgrid))::Float64
	xyidx = setdiff(1:3, zidx)
	xyidx = SVector{2, Int}(xyidx[1], xyidx[2])
	rlattice = SMatrix{2, 2}(rlattice[xyidx, xyidx])
	kgrid = kgrid[xyidx]
	return v(Val(:head), rlattice, kgrid)
end
function (v::ReciprocalCoulombKeldyshn0)(::Val{:head}, rlattice::AbstractMatrix, kgrid::Union{AbstractVector, Tuple})::Float64
	all(signbit, kgrid) && error("nk should be positive")

	minirlattice = hcat(rlattice[1:2, 1] ./ kgrid[1], rlattice[1:2, 2] ./ kgrid[2])
	miniFBZ = WignerSeitz(minirlattice)

	qsz2 = minimum(miniFBZ.offsets) / 2
	qsz = sqrt(qsz2)

	a = Vector{Float64}(undef, 2)
	b = Vector{Float64}(undef, 2)
	for i in 1:2
		(a[i], b[i]) = extrema(k -> k[i], miniFBZ.normals)
	end

	function integrand(k)
		k² = k[1] * k[1] + k[2] * k[2]
		k² < qsz2 && return 0.0
		k ∈ miniFBZ && begin
			k = sqrt(k²)
			return 1 / (1 + v.r₀ * k) / k
		end
		return 0.0
	end

	(result, err) = hcubature(integrand, a, b; norm = abs, rtol = 1e-8, atol = 1e-10)
	result_circle = 2π * log(1 + v.r₀ * qsz) / v.r₀

	nk = prod(kgrid)
	NS = nk * v.S
	return v.coff * (result + result_circle) * NS / (2π)^2
end
struct ReciprocalCoulombKeldysh0 <: ReciprocalCoulombKeldysh
	ϵ₁::Float64
	ϵ₂::Float64
	S::Float64
	coff::Float64
end
function (v::ReciprocalCoulombKeldysh0)(k::AbstractVector{<:Real}; atol::Real = 1e-10)::Float64
	k = sqrt(k ⋅ k)
	return k ≤ atol ? 0.0 : v.coff / k
end
function (v::ReciprocalCoulombKeldysh0)(k::Real; atol::Real = 1e-10)::Float64
	signbit(k) && error("k should be positive")
	return k ≤ atol ? 0.0 : v.coff / k
end
function (v::ReciprocalCoulombKeldysh0)(::Val{:k²}, k²; atol::Real = 1e-10)::Float64
	signbit(k²) && error("k² should be positive")
	return k² ≤ atol ? 0.0 : begin
		k = sqrt(k²)
		v.coff / k
	end
end
function (v::ReciprocalCoulombKeldysh0)(::Val{:kscale}, k::AbstractVector{<:Real})::Float64
	return v.coff
end
function (v::ReciprocalCoulombKeldysh0)(::Val{:kscale}, k::Real)::Float64
	signbit(k) && error("k should be positive")
	return v.coff
end
function (v::ReciprocalCoulombKeldysh0)(::Val{:head}, nk::Integer)::Float64
	signbit(nk) && error("nk should be positive")
	NS = nk * v.S
	qsz = sqrt(4π / NS)
	return v.coff * NS * qsz / 2π
end
function (v::ReciprocalCoulombKeldysh0)(::Val{:head}, rlattice::ReciprocalLattice, kgrid::Union{AbstractVector, Tuple},
	zidx::Integer = findfirst(nk->nk == 1, kgrid))::Float64
	xyidx = setdiff(1:3, zidx)
	xyidx = SVector{2, Int}(xyidx[1], xyidx[2])
	rlattice = SMatrix{2, 2}(rlattice[xyidx, xyidx])
	kgrid = kgrid[xyidx]
	return v(Val(:head), rlattice, kgrid)
end
function (v::ReciprocalCoulombKeldysh0)(::Val{:head}, rlattice::AbstractMatrix, kgrid::Union{AbstractVector, Tuple})::Float64
	all(signbit, kgrid) && error("nk should be positive")

	minirlattice = hcat(rlattice[1:2, 1] ./ kgrid[1], rlattice[1:2, 2] ./ kgrid[2])
	miniFBZ = WignerSeitz(minirlattice)

	qsz2 = minimum(miniFBZ.offsets) / 2
	qsz = sqrt(qsz2)

	a = Vector{Float64}(undef, 2)
	b = Vector{Float64}(undef, 2)
	for i in 1:2
		(a[i], b[i]) = extrema(k -> k[i], miniFBZ.normals)
	end

	function integrand(k)
		k² = k[1] * k[1] + k[2] * k[2]
		k² < qsz2 && return 0.0
		k ∈ miniFBZ && begin
			k = sqrt(k²)
			return 1 / k
		end
		return 0.0
	end

	(result, err) = hcubature(integrand, a, b; norm = abs, rtol = 1e-8, atol = 1e-10)
	result_circle = 2π * qsz

	nk = prod(kgrid)
	NS = nk * v.S
	return v.coff * (result + result_circle) * NS / (2π)^2
end
"""
	ReciprocalCoulomb(:keldysh; r₀::Real = 1, S::Real = 1, ϵ₁::Real = 1, ϵ₂::Real = 1) -> ReciprocalCoulombKeldysh

Calculate Coulomb potential generated by a point charge in 2D in reciprocal space:

```math
V(k) = \\frac{e²}{2Aϵ₀ϵdS} \\frac{1}{k(1+r₀k)}
ϵd = (ϵ₁ + ϵ₂) / 2
```

unit: k(1/Å), V(eV), r₀(Å), S(Å²), ϵ₁(dimensionless), ϵ₂(dimensionless).

# Example

```julia
julia> V = ReciprocalCoulomb(:keldysh; r₀ = 1, S = 10)
julia> k = [1.0, 0.0, 0.0] # or k = 1.0
julia> V(k)
```

We provide a method to calculate the head value:

```math
V(head) = \\frac{NS}{(2π)²} \\int_{0}^{q_{sz}} V(k) 2πk dk , 
πq_{sz}^2 = \\frac{(2π)²}{NS}
```

```julia
julia> V(:head, N)
```
"""
function ReciprocalCoulomb(::Val{:keldysh}; r₀::Real = 1, ϵ₁::Real = 1, ϵ₂::Real = 1, S::Real = 1, atol::Real = 1e-10)
	@assert r₀ ≥ 0 "r₀ should be positive"
	@assert ϵ₁ > 0 "ϵ₁ should be positive"
	@assert ϵ₂ > 0 "ϵ₂ should be positive"
	@assert S > 0 "S should be positive"
	ϵd = (ϵ₁ + ϵ₂) / 2
	coff = CoulombScale * 2π / (S * ϵd)
	if r₀ ≤ atol
		return ReciprocalCoulombKeldysh0(ϵ₁, ϵ₂, S, coff)
	else
		return ReciprocalCoulombKeldyshn0(r₀, ϵ₁, ϵ₂, S, coff)
	end
end
