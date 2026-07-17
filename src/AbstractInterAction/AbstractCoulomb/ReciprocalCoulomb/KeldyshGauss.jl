abstract type ReciprocalCoulombKeldyshGauss <: ReciprocalCoulomb end
struct ReciprocalCoulombKeldyshGaussn0 <: ReciprocalCoulombKeldyshGauss
	r₀::Float64
	ϵ₁::Float64
	ϵ₂::Float64
	α::Float64
	αz::Float64
	S::Float64
	coff::Float64
	α²4::Float64
	αz2::Float64
end
function (v::ReciprocalCoulombKeldyshGaussn0)(k::AbstractVector{<:Real}; atol::Real = 1e-10)::Float64
	k² = k ⋅ k
	k = sqrt(k²)
	return k ≤ atol ? 0.0 : v.coff * exp(-k² / v.α²4) * erfcx(k / v.αz2) / (1 + v.r₀ * k) / k
end
function (v::ReciprocalCoulombKeldyshGaussn0)(k::Real; atol::Real = 1e-10)::Float64
	signbit(k) && error("k should be positive")
	return k ≤ atol ? 0.0 : v.coff * exp(-k * k / v.α²4) * erfcx(k / v.αz2) / (1 + v.r₀ * k) / k
end
function (v::ReciprocalCoulombKeldyshGaussn0)(::Val{:k²}, k²; atol::Real = 1e-10)::Float64
	signbit(k²) && error("k² should be positive")
	return k² ≤ atol ? 0.0 : begin
		k = sqrt(k²)
		v.coff * exp(-k² / v.α²4) * erfcx(k / v.αz2) / (1 + v.r₀ * k) / k
	end
end
function (v::ReciprocalCoulombKeldyshGaussn0)(::Val{:kscale}, k::AbstractVector{<:Real})::Float64
	k² = k ⋅ k
	k = sqrt(k²)
	return v.coff * exp(-k² / v.α²4) * erfcx(k / v.αz2) / (1 + v.r₀ * k)
end
function (v::ReciprocalCoulombKeldyshGaussn0)(::Val{:kscale}, k::Real)::Float64
	signbit(k) && error("k should be positive")
	return v.coff * exp(-k * k / v.α²4) * erfcx(k / v.αz2) / (1 + v.r₀ * k)
end
function (v::ReciprocalCoulombKeldyshGaussn0)(::Val{:head}, nk::Integer)::Float64
	signbit(nk) && error("nk should be positive")
	NS = nk * v.S
	qsz = sqrt(4 * π / NS)

	integrand(k) = exp(-k * k / v.α²4) * erfcx(k / v.αz2) / (1 + v.r₀ * k)
	result, err = quadgk(integrand, 0.0, qsz; rtol = 1e-8, atol = 1e-10)

	return v.coff * result * NS / 2π
end
function (v::ReciprocalCoulombKeldyshGaussn0)(::Val{:head}, rlattice::ReciprocalLattice, kgrid::Union{AbstractVector, Tuple},
	zidx::Integer = findfirst(nk->nk == 1, kgrid))::Float64
	xyidx = setdiff(1:3, zidx)
	xyidx = SVector{2, Int}(xyidx[1], xyidx[2])
	rlattice = SMatrix{2, 2}(rlattice[xyidx, xyidx])
	kgrid = kgrid[xyidx]
	return v(Val(:head), rlattice, kgrid)
end
function (v::ReciprocalCoulombKeldyshGaussn0)(::Val{:head}, rlattice::AbstractMatrix, kgrid::Union{AbstractVector, Tuple})::Float64
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
			return exp(-k² / v.α²4) * erfcx(k / v.αz2) / (1 + v.r₀ * k) / k
		end
		return 0.0
	end

	(result, err) = hcubature(integrand, a, b; norm = abs, rtol = 1e-8, atol = 1e-10)

	integrand_circle(k) = exp(-k * k / v.α²4) * erfcx(k / v.αz2) / (1 + v.r₀ * k)
	result_circle, err = quadgk(integrand_circle, 0.0, qsz; rtol = 1e-8, atol = 1e-10)
	result_circle = 2π * result_circle

	nk = prod(kgrid)
	NS = nk * v.S
	return v.coff * (result + result_circle) * NS / (2π)^2
end
struct ReciprocalCoulombKeldyshGauss0 <: ReciprocalCoulombKeldyshGauss
	ϵ₁::Float64
	ϵ₂::Float64
	α::Float64
	αz::Float64
	S::Float64
	coff::Float64
	α²4::Float64
	αz2::Float64
end
function (v::ReciprocalCoulombKeldyshGauss0)(k::AbstractVector{<:Real}; atol::Real = 1e-10)::Float64
	k² = k ⋅ k
	k = sqrt(k²)
	return k ≤ atol ? 0.0 : v.coff * exp(-k² / v.α²4) * erfcx(k / v.αz2) / k
end
function (v::ReciprocalCoulombKeldyshGauss0)(k::Real; atol::Real = 1e-10)::Float64
	signbit(k) && error("k should be positive")
	return k ≤ atol ? 0.0 : v.coff * exp(-k * k / v.α²4) * erfcx(k / v.αz2) / k
end
function (v::ReciprocalCoulombKeldyshGauss0)(::Val{:k²}, k²; atol::Real = 1e-10)::Float64
	signbit(k²) && error("k² should be positive")
	return k² ≤ atol ? 0.0 : begin
		k = sqrt(k²)
		v.coff * exp(-k² / v.α²4) * erfcx(k / v.αz2) / k
	end
end
function (v::ReciprocalCoulombKeldyshGauss0)(::Val{:kscale}, k::AbstractVector{<:Real})::Float64
	k² = k ⋅ k
	k = sqrt(k²)
	return v.coff * exp(-k² / v.α²4) * erfcx(k / v.αz2)
end
function (v::ReciprocalCoulombKeldyshGauss0)(::Val{:kscale}, k::Real)::Float64
	signbit(k) && error("k should be positive")
	return v.coff * exp(-k * k / v.α²4) * erfcx(k / v.αz2)
end
function (v::ReciprocalCoulombKeldyshGauss0)(::Val{:head}, nk::Integer)::Float64
	signbit(nk) && error("nk should be positive")
	NS = nk * v.S
	qsz = sqrt(4π / NS)

	integrand(k) = exp(-k * k / v.α²4) * erfcx(k / v.αz2)
	result, err = quadgk(integrand, 0.0, qsz; rtol = 1e-8, atol = 1e-10)

	return v.coff * result * NS / 2π
end
function (v::ReciprocalCoulombKeldyshGauss0)(::Val{:head}, rlattice::ReciprocalLattice, kgrid::Union{AbstractVector, Tuple},
	zidx::Integer = findfirst(nk->nk == 1, kgrid))::Float64
	xyidx = setdiff(1:3, zidx)
	xyidx = SVector{2, Int}(xyidx[1], xyidx[2])
	rlattice = SMatrix{2, 2}(rlattice[xyidx, xyidx])
	kgrid = kgrid[xyidx]
	return v(Val(:head), rlattice, kgrid)
end
function (v::ReciprocalCoulombKeldyshGauss0)(::Val{:head}, rlattice::AbstractMatrix, kgrid::Union{AbstractVector, Tuple})::Float64
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
			return exp(-k * k / v.α²4) * erfcx(k / v.αz2) / k
		end
		return 0.0
	end

	(result, err) = hcubature(integrand, a, b; norm = abs, rtol = 1e-8, atol = 1e-10)

	integrand_circle(k) = exp(-k * k / v.α²4) * erfcx(k / v.αz2)
	result_circle, err = quadgk(integrand_circle, 0.0, qsz; rtol = 1e-8, atol = 1e-10)
	result_circle = 2π * result_circle

	nk = prod(kgrid)
	NS = nk * v.S
	return v.coff * (result + result_circle) * NS / (2π)^2
end
"""
	ReciprocalCoulomb(:keldysh_gauss; r₀::Real = 0, ϵ₁::Real = 1, ϵ₂::Real = 1, α::Real = 1.0, αz::Real = 1.0, S::Real = 1.0, atol::Real = 0) -> ReciprocalCoulombKeldyshGauss

Under 2D cases, calculate Coulomb potential generated by a charge of Gaussian distribution in reciprocal space:

```math
V(k) = \\frac{e²}{2ϵ₀ϵdS} \\frac{e^{-\\frac{k^2}{4α^2}} erfcx(\\frac{k}{2αz}) }{k(1+r₀k)}
ϵd = (ϵ₁ + ϵ₂) / 2
```

unit: k(1/Å), V(eV), r₀(Å), ϵ₁(dimensionless), ϵ₂(dimensionless), α(1/Å), αz(1/Å), S(Å²).

# Example

```julia
julia> V = ReciprocalCoulomb(:keldysh_gauss; r₀ = 1)
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
function ReciprocalCoulomb(::Val{:keldysh_gauss}; r₀::Real = 0, ϵ₁::Real = 1, ϵ₂::Real = 1, α::Real = 1, αz::Real = 1, S::Real = 1, atol::Real = 0)
	@assert r₀ ≥ 0 "r₀ should be positive"
	@assert ϵ₁ > 0 "ϵ₁ should be positive"
	@assert ϵ₂ > 0 "ϵ₂ should be positive"
	@assert α > 0 "α should be positive"
	@assert αz > 0 "αz should be positive"
	@assert S > 0 "S should be positive"
	ϵd = (ϵ₁ + ϵ₂) / 2
	coff = CoulombScale * 2π / (S * ϵd)
	α²4 = 4 * α * α
	αz2 = 2 * αz
	if r₀ ≤ atol
		return ReciprocalCoulombKeldyshGauss0(ϵ₁, ϵ₂, α, αz, S, coff, α²4, αz2)
	else
		return ReciprocalCoulombKeldyshGaussn0(r₀, ϵ₁, ϵ₂, α, αz, S, coff, α²4, αz2)
	end
end
