struct LongRangePartKeldyshGauss <: LongRangePart
	# wannier structure.
	lattice::Lattice{Float64}
	rlattice::ReciprocalLattice{Float64}
	period::Vec3{Bool}
	xyindex::SVector{2, Int}
	zindex::Int
	lattice2D::SMatrix{2, 2, Float64, 4}
	rlattice2D::SMatrix{2, 2, Float64, 4}
	rldot::SMatrix{2, 2, Float64, 4}
	S::Float64
	norb::Int
	orblocat::Vector{SVector{2, Float64}}
	# Coulomb information.
	r₀::Float64
	ϵ₁::Float64
	ϵ₂::Float64
	α::Matrix{Float64}
	αz::Matrix{Float64}
	φ::Matrix{ReciprocalCoulombKeldyshGauss}
	coff::Float64
	ϵd::Float64
	α²4::Matrix{Float64}
	αz2::Matrix{Float64}
	# G grid for the summation of long-range correction.
	nG::Int
	# only contain nonzero G vectors, the length equal to nG-1.
	nGminus1::Int
	G_frac::Vector{SVector{2, Int}}
	G_car::Vector{SVector{2, Float64}}
	G_norm2::Vector{Float64}
	Δorb::Matrix{SVector{2, Float64}}
	Δorb_car::Matrix{SVector{2, Float64}}
	G_orb_phase::Array{ComplexF64, 3}
	Δorb_norm::Matrix{Float64}
end
@inline (v::LongRangePartKeldyshGauss)(::Val{+}, A, k::ReducedCoordinates, args...; kwargs...) =
	v(Val(+), A, k[v.xyindex], args...; kwargs...)
function (v::LongRangePartKeldyshGauss)(::Val{+}, A, k::SVector{2, <:Real}; isΓ::Bool = false, nearΓ::Bool = isΓ)
	if nearΓ
		return _LongRangePartKeldyshGauss_wok!(A, v, normalize_kdirect(k), Val(isΓ))
	else
		return _LongRangePartKeldyshGauss_wik!(A, v, normalize_kdirect(k))
	end
end
function (v::LongRangePartKeldyshGauss)(::Val{+}, A, k::SVector{2, <:Real}, nkorkgrid; isΓ::Bool = false, nearΓ::Bool = isΓ)
	if nearΓ
		_LongRangePartKeldyshGauss_head!(A, v, nkorkgrid)
		return _LongRangePartKeldyshGauss_wok!(A, v, normalize_kdirect(k), Val(isΓ))
	else
		return _LongRangePartKeldyshGauss_wik!(A, v, normalize_kdirect(k))
	end
end
(v::LongRangePartKeldyshGauss)(::Val{:G0scale}, A, k::ReducedCoordinates) = v(Val(:G0scale), A, k[v.xyindex])
function (v::LongRangePartKeldyshGauss)(::Val{:G0scale}, A, k::SVector{2, <:Real})
	for j in 2:v.norb, i in 1:(j-1)
		A[i, j] = v.φ[i, j](Val(:kscale), k) * cispi(2 * (k ⋅ v.Δorb[i, j]))
		A[j, i] = conj(A[i, j])
	end
	for i in 1:v.norb
		A[i, i] = v.φ[i, i](Val(:kscale), k)
	end
	return A
end
# function (v::LongRangePartKeldyshGauss)(::Val{:taylor_n1})
# 	return v.coff
# end
function (v::LongRangePartKeldyshGauss)(::Val{:taylor_0})
	return - v.coff .* (v.r₀ .+ 1.0 ./ (√π .* v.αz))
end
function (v::LongRangePartKeldyshGauss)(::Val{:taylor_1})
	return v.coff .* (v.r₀^2 .- 1.0 ./ v.α²4 .+ 1.0 ./ (v.αz2 .* v.αz2) .+ r₀ ./ (√π .* v.αz))
end
function _LongRangePartKeldyshGauss_wik!(A, v, k)

	k_car = v.rlattice2D * k
	kG_norm2 = map(v.G_car) do G
		kG = k_car + G
		kG ⋅ kG
	end
	k_norm2 = k_car ⋅ k_car

	φkG = Vector{Float64}(undef, v.nGminus1)
	for j in 2:v.norb, i in 1:(j-1)
		φᵢⱼ = v.φ[i, j]
		map(enumerate(kG_norm2)) do (iG, kG2)
			φkG[iG] = φᵢⱼ(Val(:k²), kG2; atol = 0.0)
		end
		vlr = cispi(2 * (k ⋅ v.Δorb[i, j])) * (sum(iG -> φkG[iG] * v.G_orb_phase[iG, i, j], Base.OneTo(v.nGminus1)) +
											   φᵢⱼ(Val(:k²), k_norm2; atol = 0.0))
		A[i, j] += vlr
		A[j, i] += conj(vlr)
	end
	for i in 1:v.norb
		φᵢᵢ = v.φ[i, i]
		map(enumerate(kG_norm2)) do (iG, kG2)
			φkG[iG] = φᵢᵢ(Val(:k²), kG2; atol = 0.0)
		end
		vlr = sum(φkG) + φᵢᵢ(Val(:k²), k_norm2; atol = 0.0)
		A[i, i] += vlr
	end

	return A
end
function _LongRangePartKeldyshGauss_wok!(A, v, k, ::Val{false})

	k_car = v.rlattice2D * k
	kG_norm2 = map(v.G_car) do G
		kG = k_car + G
		kG ⋅ kG
	end

	φkG = Vector{Float64}(undef, v.nGminus1)
	for j in 2:v.norb, i in 1:(j-1)
		φᵢⱼ = v.φ[i, j]
		map(enumerate(kG_norm2)) do (iG, kG2)
			φkG[iG] = φᵢⱼ(Val(:k²), kG2; atol = 0.0)
		end
		vlr = cispi(2 * (k ⋅ v.Δorb[i, j])) * sum(iG -> φkG[iG] * v.G_orb_phase[iG, i, j], Base.OneTo(v.nGminus1))
		A[i, j] += vlr
		A[j, i] += conj(vlr)
	end
	for i in 1:v.norb
		φᵢᵢ = v.φ[i, i]
		map(enumerate(kG_norm2)) do (iG, kG2)
			φkG[iG] = φᵢᵢ(Val(:k²), kG2; atol = 0.0)
		end
		vlr = sum(φkG)
		A[i, i] += vlr
	end

	return A
end
function _LongRangePartKeldyshGauss_wok!(A, v, k, ::Val{true})

	φG = Vector{Float64}(undef, v.nGminus1)
	for j in 2:v.norb, i in 1:(j-1)
		φᵢⱼ = v.φ[i, j]
		map(enumerate(v.G_norm2)) do (iG, G2)
			φG[iG] = φᵢⱼ(Val(:k²), G2; atol = 0.0)
		end
		vlr = sum(iG -> φG[iG] * real(v.G_orb_phase[iG, i, j]), Base.OneTo(v.nGminus1))
		A[i, j] += vlr
		A[j, i] += vlr
	end
	for i in 1:v.norb
		φᵢᵢ = v.φ[i, i]
		map(enumerate(v.G_norm2)) do (iG, G2)
			φG[iG] = φᵢᵢ(Val(:k²), G2; atol = 0.0)
		end
		vlr = sum(φG)
		A[i, i] += vlr
	end

	return A
end
function _LongRangePartKeldyshGauss_head!(A, v, nk::Integer)
	signbit(nk) && error("nk should be positive")
	NS = nk * v.S
	qsz = sqrt(4 * π / NS)
	coff = CoulombScale * nk / (2π * v.ϵd)

	for j in 2:v.norb, i in 1:(j-1)
		integrand(k) = exp(-k * k / v.α²4[i, j]) * erfcx(k / v.αz2[i, j]) / (1 + v.r₀ * k)
		result, err = quadgk(integrand, 0.0, qsz; rtol = 1e-8, atol = 1e-10)
		head = coff * 2π * result
		A[i, j] += head
		A[j, i] += head
	end
	for i in 1:v.norb
		integrand(k) = exp(-k * k / v.α²4[i, i]) * erfcx(k / v.αz2[i, i]) / (1 + v.r₀ * k)
		result, err = quadgk(integrand, 0.0, qsz; rtol = 1e-8, atol = 1e-10)
		head = coff * 2π * result
		A[i, i] += head
	end

	return A
end
function _LongRangePartKeldyshGauss_head!(A, v, kgrid::Union{AbstractVector{<:Integer}, NTuple{3, <:Integer}})
	all(signbit, kgrid) && error("nk should be positive")
	length(kgrid) == 3 || error("kgrid should have three components")
	kgrid[v.zindex] == 1 || error("kgrid should have only one point along the non-periodic direction")
	kgrid = kgrid[v.xyindex]

	minirlattice = [v.rlattice2D[1:2, 1] ./ kgrid[1] v.rlattice2D[1:2, 2] ./ kgrid[2]]
	miniFBZ = WignerSeitz(minirlattice)

	qsz2 = minimum(miniFBZ.offsets) / 2
	qsz = sqrt(qsz2)

	a = Vector{Float64}(undef, 2)
	b = Vector{Float64}(undef, 2)
	for i in 1:2
		(a[i], b[i]) = extrema(k -> k[i], miniFBZ.normals)
	end

	nk = prod(kgrid)
	coff = CoulombScale * nk / (2π * v.ϵd)

	for j in 2:v.norb, i in 1:(j-1)
		function integrand(k)
			k² = k[1] * k[1] + k[2] * k[2]
			k² < qsz2 && return 0.0
			k ∈ miniFBZ && begin
				k = sqrt(k²)
				return exp(-k² / v.α²4[i, j]) * erfcx(k / v.αz2[i, j]) / (1 + v.r₀ * k) / k
			end
			return 0.0
		end
		(result, err) = hcubature(integrand, a, b; norm = abs, rtol = 1e-8, atol = 1e-10)
		integrand_circle(k) = exp(-k * k / v.α²4[i, j]) * erfcx(k / v.αz2[i, j]) / (1 + v.r₀ * k)
		result_circle, err = quadgk(integrand_circle, 0.0, qsz; rtol = 1e-8, atol = 1e-10)
		head = coff * (2π * result_circle + result)
		A[i, j] += head
		A[j, i] += head
	end
	for i in 1:v.norb
		function integrand(k)
			k² = k[1] * k[1] + k[2] * k[2]
			k² < qsz2 && return 0.0
			k ∈ miniFBZ && begin
				k = sqrt(k²)
				return exp(-k² / v.α²4[i, i]) * erfcx(k / v.αz2[i, i]) / (1 + v.r₀ * k) / k
			end
			return 0.0
		end
		(result, err) = hcubature(integrand, a, b; norm = abs, rtol = 1e-8, atol = 1e-10)
		integrand_circle(k) = exp(-k * k / v.α²4[i, i]) * erfcx(k / v.αz2[i, i]) / (1 + v.r₀ * k)
		result_circle, err = quadgk(integrand_circle, 0.0, qsz; rtol = 1e-8, atol = 1e-10)
		head = coff * (2π * result_circle + result)
		A[i, i] += head
	end

	return A
end
function LongRangePart(::Val{:keldysh_gauss}, lattice::Lattice,
	orblocat::Vector{<:ReducedCoordinates}; norb::Integer = length(orblocat),
	period = Bool[1, 1, 0],
	α::Real = 1, αmat::AbstractMatrix{<:Real} = fill(α, (norb, norb)),
	αz::Real = 1, αzmat::AbstractMatrix{<:Real} = fill(αz, (norb, norb)),
	r₀::Real = 0, ϵ₁::Real = 1, ϵ₂::Real = 1, r₀atol::Real = 0, δ::Real = 1e-8)

	period = Vec3{Bool}(period)
	@assert count(period) == 2 "Exactly two dimensions should be periodic."
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

	@assert all(iszero, lattice[xyindex, zindex]) "Wrong lattice: the non-periodic direction should be perpendicular to the periodic plane."
	@assert all(iszero, lattice[zindex, xyindex]) "Wrong lattice: the non-periodic direction should be perpendicular to the periodic plane."

	rlattice = reciprocal(lattice)
	lattice2D = SMatrix{2, 2, Float64, 4}(lattice[xyindex, xyindex])
	rlattice2D = SMatrix{2, 2, Float64, 4}(rlattice[xyindex, xyindex])
	rldot = SMatrix{2, 2, Float64, 4}(transpose(rlattice2D) * rlattice2D)

	𝐚 = lattice2D[:, 1]
	𝐛 = lattice2D[:, 2]
	S = abs(𝐚[1] * 𝐛[2] - 𝐚[2] * 𝐛[1])

	orblocat = map(x->SVector{2, Float64}(x[xyindex]), orblocat)

	αmat, αzmat = _LongRangePartKeldyshGauss_symmetrize_αmat(αmat, αzmat)
	φ = [ReciprocalCoulomb(:keldysh_gauss; r₀, ϵ₁, ϵ₂, α = αmat[i, j], αz = αzmat[i, j], S, atol = r₀atol) for i in 1:norb, j in 1:norb]
	ϵd = (ϵ₁ + ϵ₂) / 2
	coff = CoulombScale * 2π / (S * ϵd)
	α²4 = 4 .* αmat .* αmat
	αz2 = 2 .* αzmat

	Gnorm_max = sqrt(-4 * maximum(αmat)^2 * log(δ)) + norm(rlattice2D)
	G_frac = gridindex(rlattice, Gnorm_max, xyindex, Val(2))
	G_frac = map(G -> SVector{2, Int}(G[xyindex]), G_frac)
	nG = length(G_frac)

	G0idx = findfirst(iszero, G_frac)
	Gn0idx = setdiff(1:nG, G0idx)
	nGminus1 = nG - 1
	G_frac = G_frac[Gn0idx]
	G_car = map(G -> rlattice2D * G, G_frac)
	G_norm2 = map(G -> G[1] * G[1] + G[2] * G[2], G_car)

	Δorb = Matrix{SVector{2, Float64}}(undef, norb, norb)
	Δorb_car = Matrix{SVector{2, Float64}}(undef, norb, norb)
	G_orb_phase = Array{ComplexF64}(undef, nGminus1, norb, norb)
	Δorb_norm = Matrix{Float64}(undef, norb, norb)
	for i in 1:norb, j in 1:i
		dorb = orblocat[i] - orblocat[j]
		dorb_car = lattice2D * dorb
		Δorb[i, j] = dorb
		Δorb[j, i] = -dorb
		Δorb_car[i, j] = dorb_car
		Δorb_car[j, i] = -dorb_car
		for (Gi, G) in enumerate(G_frac)
			G_orb_phase[Gi, i, j] = cispi(2 * (G ⋅ Δorb[i, j]))
		end
		G_orb_phase[:, j, i] .= conj.(G_orb_phase[:, i, j])
		dorb_norm = norm(dorb_car)
		Δorb_norm[i, j] = dorb_norm
		Δorb_norm[j, i] = Δorb_norm[i, j]
	end

	return LongRangePartKeldyshGauss(lattice, rlattice, period, xyindex, zindex,
		lattice2D, rlattice2D, rldot, S, norb, orblocat,
		r₀, ϵ₁, ϵ₂, αmat, αzmat, φ, coff, ϵd, α²4, αz2,
		nG, nGminus1, G_frac, G_car, G_norm2,
		Δorb, Δorb_car, G_orb_phase, Δorb_norm)
end
function _LongRangePartKeldyshGauss_symmetrize_αmat(αmat, αzmat)
	norb = size(αmat, 1)
	αmat_sym = Matrix{Float64}(undef, norb, norb)
	αzmat_sym = Matrix{Float64}(undef, norb, norb)
	for j in 2:norb, i in 1:(j-1)
		αmat_sym[i, j] = αmat[i, j]
		αmat_sym[j, i] = αmat[i, j]
		αzmat_sym[i, j] = αzmat[i, j]
		αzmat_sym[j, i] = αzmat[i, j]
	end
	for i in 1:norb
		αmat_sym[i, i] = αmat[i, i]
		αzmat_sym[i, i] = αzmat[i, i]
	end
	return αmat_sym, αzmat_sym
end
