#Note ϕK(0) is redefined to 0, is not divergent.
#But we need its divergence, so we can't input k=0.
export GaussLRCorrection, GaussLRCorrection3D, GaussLRCorrection2D, GaussLRCorrection1D, GaussLRCorrection0D

abstract type GaussLRCorrection <: AbstractLRCorrection end

function (v::GaussLRCorrection)(k::AbstractVector)
	A = zeros(ComplexF64, v.norb, v.norb)
	return v(Val(+), A, k)
end
function (v::GaussLRCorrection)(A::AbstractMatrix, k::AbstractVector)
	size(A) == (v.norb, v.norb) || error("Buffer size mismatch.")
	A .= 0
	return v(Val(+), A, k)
end
function (v::GaussLRCorrection)(k::AbstractVector, nk::Integer)
	A = zeros(ComplexF64, v.norb, v.norb)
	return v(Val(+), A, k, nk)
end
function (v::GaussLRCorrection)(A::AbstractMatrix, k::AbstractVector, nk::Integer)
	size(A) == (v.norb, v.norb) || error("Buffer size mismatch.")
	A .= 0
	return v(Val(+), A, k, nk)
end

struct GaussLRCorrection3D <: GaussLRCorrection
	lattice::Lattice{Float64}
	ϵ::Float64
	α::Float64
	Ω::Float64
	α²4::Float64
	φ::ReciprocalGauss
	rlattice::ReciprocalLattice{Float64}
	NG::Int
	G_frac::Vector{ReducedCoordinates{Int}}
	G_car::Vector{CartesianCoordinates{Float64}}
	G0idx::Int
	Gn0idx::Vector{Int}
	norb::Int
	Δorb::Matrix{ReducedCoordinates{Float64}}
	G_orb_phase::Array{ComplexF64, 3}
	Δorb_norm::Matrix{Float64}
	Δorb_norm_divπ::Matrix{Float64}
end
function (v::GaussLRCorrection3D)(::Val{+}, A, k::AbstractVector, nk::Integer)
	if iszero(k)
		return v(Val(+), A, Val(:head), nk)
	else
		return v(Val(+), A, k)
	end
end
function (v::GaussLRCorrection3D)(::Val{+}, A, k::AbstractVector)

	k_car = v.rlattice * ReducedCoordinates(k)
	φkG = map(G -> v.φ(k_car + G), v.G_car)

	for j in 2:v.norb, i in 1:j-1
		vlr = cispi(2 * (k ⋅ v.Δorb[i, j])) * sum(Gi -> φkG[Gi] * v.G_orb_phase[Gi, i, j], Base.OneTo(v.NG))
		A[i, j] += vlr
		A[j, i] += conj(vlr)
	end
	Vᵢᵢ = sum(φkG)
	for i in 1:v.norb
		A[i, i] += Vᵢᵢ
	end

	return A
end
function (v::GaussLRCorrection3D)(::Val{+}, A, k::ReducedCoordinates)

	k_car = v.rlattice * k
	φkG = map(G -> v.φ(k_car + G), v.G_car)

	for j in 2:v.norb, i in 1:j-1
		vlr = cispi(2 * (k ⋅ v.Δorb[i, j])) * sum(Gi -> φkG[Gi] * v.G_orb_phase[Gi, i, j], Base.OneTo(v.NG))
		A[i, j] += vlr
		A[j, i] += conj(vlr)
	end
	Vᵢᵢ = sum(φkG)
	for i in 1:v.norb
		A[i, i] += Vᵢᵢ
	end

	return A
end
function (v::GaussLRCorrection3D)(::Val{+}, A, ::Val{:head}, nk::Integer)
	signbit(nk) && error("nk should be positive")
	NΩ = nk * v.Ω
	qsz = (6 * π^2 / NΩ)^(1 // 3)
	T = CoulombScale * 2 * nk / π / v.ϵ

	φG = map(G -> v.φ(G), v.G_car)

	for j in 2:v.norb, i in 1:j-1
		vlr = real(sum(Gi -> φG[Gi] * v.G_orb_phase[Gi, i, j], v.Gn0idx))
		integrand(k) = exp(-k * k / v.α²4) * sinc(k * v.Δorb_norm_divπ[i, j])
		result, err = quadgk(integrand, 0.0, qsz, rtol = 1e-9)
		head = result * T
		A[i, j] += vlr + head
		A[j, i] += vlr + head
	end
	Vᵢᵢ = sum(φG)
	φhead = v.φ(Val(:head), nk)
	for i in 1:v.norb
		A[i, i] += Vᵢᵢ + φhead
	end

	return A
end
function GaussLRCorrection3D(lattice::Lattice, orblocat_frac::AbstractVector, α::Real; δ = 1e-8, ϵ = 1)

	G2_max = -4 * α^2 * log(δ) * 1.2

	rlattice = reciprocal(lattice)
	b₁ = rlattice[:, 1]
	b₂ = rlattice[:, 2]
	b₃ = rlattice[:, 3]

	V_BZ = ((b₁ × b₂) ⋅ b₃)

	h₁ = V_BZ / norm(b₂ × b₃)
	h₂ = V_BZ / norm(b₃ × b₁)
	h₃ = V_BZ / norm(b₁ × b₂)

	Ggrid = Int.(cld.(√G2_max * 1.1, [h₁, h₂, h₃])) * 2 .+ 1
	Ggrid = gridindex(Ggrid)

	G_car = map(G -> rlattice * G, Ggrid)
	Gnorm2 = map(G -> G[1] * G[1] + G[2] * G[2] + G[3] * G[3], G_car)

	I = Gnorm2 .< G2_max
	G_frac = Ggrid[I]
	G_car = G_car[I]
	NG = length(G_frac)

	Ω = abs((lattice[:, 1] × lattice[:, 2]) ⋅ lattice[:, 3])
	φ = ReciprocalGauss(; ϵ, α, Ω)

	norb = length(orblocat_frac)
	Δorb = Matrix{ReducedCoordinates{Float64}}(undef, norb, norb)
	G_orb_phase = Array{ComplexF64}(undef, NG, norb, norb)
	Δorb_norm = Matrix{Float64}(undef, norb, norb)
	Δorb_norm_divπ = Matrix{Float64}(undef, norb, norb)
	for i in 1:norb, j in 1:i
		dorb = orblocat_frac[i] - orblocat_frac[j]
		Δorb[i, j] = dorb
		Δorb[j, i] = -dorb
		for (Gi, G) in enumerate(G_frac)
			G_orb_phase[Gi, i, j] = cispi(2 * (G ⋅ Δorb[i, j]))
		end
		G_orb_phase[:, j, i] .= conj.(G_orb_phase[:, i, j])
		Δorb_norm[i, j] = norm(lattice * dorb)
		Δorb_norm[j, i] = Δorb_norm[i, j]
		Δorb_norm_divπ[i, j] = Δorb_norm[i, j] / π
		Δorb_norm_divπ[j, i] = Δorb_norm_divπ[i, j]
	end

	G0idx = findfirst(iszero, G_frac)
	Gn0idx = setdiff(1:NG, G0idx)

	α²4 = 4 * α^2

	return GaussLRCorrection3D(lattice, ϵ, α, Ω, α²4, φ, rlattice,
		NG, G_frac, G_car, G0idx, Gn0idx,
		norb, Δorb, G_orb_phase, Δorb_norm, Δorb_norm_divπ)
end
struct GaussLRCorrection2D <: GaussLRCorrection
	φ::ReciprocalGauss2D
	rlattice::ReciprocalLattice
	NG::Int
	G_frac::Vector{Vec3{Int}}
	G_car::Vector{Vec3{Float64}}
	G_orb_phase::Array{ComplexF64}
	norb::Int
	Δorb::Matrix{Vec3{Float64}}
	Δorbz::Matrix{Float64}
	xyindex::SVector{2, Int}
	zindex::Int
end
function (v::GaussLRCorrection2D)(::Val{+}, A, k::AbstractVector, nk::Integer)
	#Only use kx,ky.
	if iszero(k[v.xyindex])
		return v(Val(+), A, Val(:head), nk)
	else
		return v(Val(+), A, k)
	end
end
function (v::GaussLRCorrection2D)(::Val{+}, A, k::AbstractVector)
	#Only use kx,ky.
	k2D = [0.0, 0.0, 0.0]
	k2D[v.xyindex] .= k[v.xyindex]

	k_car = v.rlattice * k2D
	kG_car = map(G -> k_car + G, v.G_car)

	for j in 2:v.norb, i in 1:j-1
		vlr = cispi(-2 * (k2D ⋅ v.Δorb[i, j])) * sum(Gi -> v.φ(kG_car[Gi], v.Δorbz[i, j]) * v.G_orb_phase[Gi, i, j], Base.OneTo(v.NG))
		A[i, j] += vlr
		A[j, i] += conj(vlr)
	end
	Vᵢᵢ = sum(Gi -> v.φ(kG_car[Gi], 0), Base.OneTo(v.NG))
	for i in 1:v.norb
		A[i, i] += Vᵢᵢ
	end

	return A
end
function (v::GaussLRCorrection2D)(::Val{+}, A, ::Val{:head}, nk::Integer)

	G0 = findfirst(iszero, v.G_frac)
	I = setdiff(1:v.NG, G0)

	for j in 2:v.norb, i in 1:j-1
		φhead = v.φ(Val(:head), nk, v.Δorbz[i, j])
		vlr = sum(Gi -> v.φ(v.G_car[Gi], v.Δorbz[i, j]) * v.G_orb_phase[Gi, i, j], I) + φhead
		A[i, j] += vlr
		A[j, i] += conj(vlr)
	end
	φhead = v.φ(Val(:head), nk, 0)
	Vᵢᵢ = sum(Gi -> v.φ(v.G_car[Gi], 0), I) + φhead
	for i in 1:v.norb
		A[i, i] += Vᵢᵢ
	end

	return A
end
function GaussLRCorrection2D(lattice::Lattice, orblocat_frac::AbstractVector, period::AbstractVector, α::Real; δ = 1e-6, ϵ = 1)

	G2_max = -4 * α^2 * log(δ) * 1.5

	rlattice = reciprocal(lattice)
	b₁ = rlattice[:, 1]
	b₂ = rlattice[:, 2]
	b₃ = rlattice[:, 3]

	V_BZ = ((b₁ × b₂) ⋅ b₃)

	h₁ = V_BZ / norm(b₂ × b₃)
	h₂ = V_BZ / norm(b₃ × b₁)
	h₃ = V_BZ / norm(b₁ × b₂)

	Ggrid = Int.(cld.(√G2_max * 1.2, [h₁, h₂, h₃])) * 2 .+ 1

	a₁ = lattice[:, 1]
	a₂ = lattice[:, 2]
	a₃ = lattice[:, 3]

	if !period[1]
		Ggrid[1] = 1
		S = norm(a₂ × a₃)
		xyindex = SVector{2, Int}(2, 3)
		zindex = 1
	elseif !period[2]
		Ggrid[2] = 1
		S = norm(a₃ × a₁)
		xyindex = SVector{2, Int}(1, 3)
		zindex = 2
	elseif !period[3]
		Ggrid[3] = 1
		S = norm(a₁ × a₂)
		xyindex = SVector{2, Int}(1, 2)
		zindex = 3
	end

	Ggrid = gridindex(Ggrid)
	G_frac = filter(G -> sum(abs2, rlattice * G) < G2_max, Ggrid)
	G_car = map(G -> rlattice * G, G_frac)

	φ = ReciprocalGauss2D(; ϵ, α, S)

	norb = length(orblocat_frac)
	Δorb = Matrix{Vec3{Float64}}(undef, norb, norb)
	NG = length(G_frac)
	G_orb_phase = Array{ComplexF64}(undef, NG, norb, norb)
	for i in 1:norb, j in 1:i
		dorb = orblocat_frac[j] - orblocat_frac[i]
		Δorb[i, j] = dorb
		Δorb[j, i] = -dorb
		for (Gi, G) in enumerate(G_frac)
			G_orb_phase[Gi, i, j] = cispi(-2 * (G ⋅ Δorb[i, j]))
		end
		G_orb_phase[:, j, i] = conj.(G_orb_phase[:, i, j])
	end

	Δorbz = map(dorb -> (lattice*dorb)[zindex], Δorb)

	return GaussLRCorrection2D(φ, rlattice, NG, G_frac, G_car, G_orb_phase, norb, Δorb, Δorbz, xyindex, zindex)
end

#TODO
struct GaussLRCorrection1D <: GaussLRCorrection
	norb::Int
end
function (v::GaussLRCorrection1D)(::Val{+}, A, para...)
	error("TODO")
	return A
end
struct GaussLRCorrection0D <: GaussLRCorrection
	norb::Int
end
function (v::GaussLRCorrection0D)(::Val{+}, A, para...)
	return A
end
