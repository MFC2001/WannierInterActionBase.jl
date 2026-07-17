struct LongRangePartGauss <: LongRangePart
	# wannier structure.
	lattice::Lattice{Float64}
	rlattice::ReciprocalLattice{Float64}
	rldot::Mat3{Float64}
	Ω::Float64
	norb::Int
	orblocat::Vector{ReducedCoordinates{Float64}}
	# Coulomb information.
	ϵ::Float64
	α::Matrix{Float64}
	φ::Matrix{ReciprocalCoulombGauss}
	coff::Float64
	α²4::Matrix{Float64}
	# G grid for the summation of long-range correction.
	nG::Int
	# only contain nonzero G vectors, the length equal to nG-1.
	nGminus1::Int
	G_frac::Vector{ReducedCoordinates{Int}}
	G_car::Vector{CartesianCoordinates{Float64}}
	G_norm2::Vector{Float64}
	Δorb::Matrix{ReducedCoordinates{Float64}}
	Δorb_car::Matrix{CartesianCoordinates{Float64}}
	G_orb_phase::Array{ComplexF64, 3}
	Δorb_norm_time2α_divπ::Matrix{Float64}
end
function (v::LongRangePartGauss)(::Val{+}, A, k::ReducedCoordinates; isΓ::Bool = false, nearΓ::Bool = isΓ)
	if nearΓ
		return _LongRangePartGauss_wok!(A, v, normalize_kdirect(k), Val(isΓ))
	else
		return _LongRangePartGauss_wik!(A, v, normalize_kdirect(k))
	end
end
function (v::LongRangePartGauss)(::Val{+}, A, k::ReducedCoordinates, nkorkgrid; isΓ::Bool = false, nearΓ::Bool = isΓ)
	if nearΓ
		_LongRangePartGauss_head!(A, v, nkorkgrid)
		return _LongRangePartGauss_wok!(A, v, normalize_kdirect(k), Val(isΓ))
	else
		return _LongRangePartGauss_wik!(A, v, normalize_kdirect(k))
	end
end
function (v::LongRangePartGauss)(::Val{:G0scale}, A, k::ReducedCoordinates)
	for j in 2:v.norb, i in 1:(j-1)
		A[i, j] = v.φ[i, j](Val(:k²scale), k) * cispi(2 * (k ⋅ v.Δorb[i, j]))
		A[j, i] = conj(A[i, j])
	end
	for i in 1:v.norb
		A[i, i] = v.φ[i, i](Val(:k²scale), k)
	end
	return A
end
# function (v::LongRangePartGauss)(::Val{:taylor_n2})
# 	return v.coff
# end
function (v::LongRangePartGauss)(::Val{:taylor_0})
	return - v.coff ./ v.α²4
end
function _LongRangePartGauss_wik!(A, v, k)

	k_car = v.rlattice * k
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
function _LongRangePartGauss_wok!(A, v, k, ::Val{false})

	k_car = v.rlattice * k
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
function _LongRangePartGauss_wok!(A, v, k, ::Val{true})

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
function _LongRangePartGauss_head!(A, v, nk::Integer)
	signbit(nk) && error("nk should be positive")

	NΩ = nk * v.Ω
	qsz = (6 * π^2 / NΩ)^(1 // 3)
	coff = CoulombScale * nk / (2 * π^2 * v.ϵ)

	for j in 2:v.norb, i in 1:(j-1)
		head = coff * 4π * √π * v.α[i, j] * erf(qsz / (2 * v.α[i, j]))
		A[i, j] += head
		A[j, i] += head
	end
	for i in 1:v.norb
		head = coff * 4π * √π * v.α[i, i] * erf(qsz / (2 * v.α[i, i]))
		A[i, i] += head
	end

	return A
end
# function _LongRangePartGauss_head!(A, v, nk::Integer)
# 	signbit(nk) && error("nk should be positive")

# 	NΩ = nk * v.Ω
# 	qsz = (6 * π^2 / NΩ)^(1 // 3)
# 	coff = CoulombScale * 2 * nk / π / v.ϵ

# 	for j in 2:v.norb, i in 1:(j-1)
# 		integrand(t) = exp(-t * t) * sinc(t * v.Δorb_norm_time2α_divπ[i, j])
# 		result, err = quadgk(integrand, 0.0, qsz / (2 * v.α[i, j]); rtol = 1e-8, atol = 1e-10)
# 		head = result * coff * 2 * v.α[i, j]
# 		A[i, j] += head
# 		A[j, i] += head
# 	end
# 	for i in 1:v.norb
# 		head = coff * √π * v.α[i, i] * erf(qsz / (2 * v.α[i, i]))
# 		A[i, i] += head
# 	end

# 	return A
# end
function _LongRangePartGauss_head!(A, v, kgrid::Union{AbstractVector{<:Integer}, NTuple{3, <:Integer}})
	all(signbit, kgrid) && error("nk should be positive")

	minirlattice = ReciprocalLattice(v.rlattice[1:3, 1] ./ kgrid[1], v.rlattice[1:3, 2] ./ kgrid[2], v.rlattice[1:3, 3] ./ kgrid[3])
	miniFBZ = WignerSeitz(minirlattice)

	qsz2 = minimum(miniFBZ.offsets) / 2
	qsz = sqrt(qsz2)

	a = Vector{Float64}(undef, 3)
	b = Vector{Float64}(undef, 3)
	for i in 1:3
		(a[i], b[i]) = extrema(k -> k[i], miniFBZ.normals)
	end

	nk = prod(kgrid)
	coff = CoulombScale * nk / (2 * π^2 * v.ϵ)

	for j in 2:v.norb, i in 1:(j-1)
		function integrand(k)
			k² = k[1] * k[1] + k[2] * k[2] + k[3] * k[3]
			k² < qsz2 && return 0.0
			k ∈ miniFBZ && return exp(-k² / v.α²4[i, j]) / k²
			return 0.0
		end
		(result, err) = hcubature(integrand, a, b; norm = abs, rtol = 1e-8, atol = 1e-10)
		head = coff * (4π * √π * v.α[i, j] * erf(qsz / (2 * v.α[i, j])) + result)
		A[i, j] += head
		A[j, i] += head
	end
	for i in 1:v.norb
		function integrand(k)
			k² = k[1] * k[1] + k[2] * k[2] + k[3] * k[3]
			k² < qsz2 && return 0.0
			k ∈ miniFBZ && return exp(-k² / v.α²4[i, i]) / k²
			return 0.0
		end
		(result, err) = hcubature(integrand, a, b; norm = abs, rtol = 1e-8, atol = 1e-10)
		head = coff * (4π * √π * v.α[i, i] * erf(qsz / (2 * v.α[i, i])) + result)
		A[i, i] += head
	end

	return A
end
# function _LongRangePartGauss_head!(A, v, kgrid::Union{AbstractVector{<:Integer}, NTuple{3, <:Integer}})
# 	all(signbit, kgrid) && error("nk should be positive")

# 	minirlattice = ReciprocalLattice(v.rlattice[1:3, 1] ./ kgrid[1], v.rlattice[1:3, 2] ./ kgrid[2], v.rlattice[1:3, 3] ./ kgrid[3])
# 	miniFBZ = WignerSeitz(minirlattice)

# 	qsz2 = minimum(miniFBZ.offsets) / 2
# 	qsz = sqrt(qsz2)

# 	a = Vector{Float64}(undef, 3)
# 	b = Vector{Float64}(undef, 3)
# 	for i in 1:3
# 		(a[i], b[i]) = extrema(k -> k[i], miniFBZ.normals)
# 	end

# 	nk = prod(kgrid)
# 	coff = CoulombScale * 2 * nk / π / v.ϵ

# 	for j in 2:v.norb, i in 1:(j-1)
# 		function integrand(k)
# 			k² = k[1] * k[1] + k[2] * k[2] + k[3] * k[3]
# 			k² < qsz2 && return 0.0
# 			k ∈ miniFBZ && return cis(k ⋅ v.Δorb_car[i, j]) * exp(-k² / v.α²4[i, j]) / k²
# 			return 0.0
# 		end
# 		(result, err) = hcubature(integrand, a, b; norm = abs, rtol = 1e-8, atol = 1e-10)
# 		integrand_sphere(t) = exp(-t * t) * sinc(t * v.Δorb_norm_time2α_divπ[i, j])
# 		result_sphere, err = quadgk(integrand_sphere, 0.0, qsz / (2 * v.α[i, j]); rtol = 1e-8, atol = 1e-10)
# 		head = coff * (2 * v.α[i, j] * result_sphere + result / (4π))
# 		A[i, j] += head
# 		A[j, i] += head
# 	end
# 	for i in 1:v.norb
# 		function integrand(k)
# 			k² = k[1] * k[1] + k[2] * k[2] + k[3] * k[3]
# 			k² < qsz2 && return 0.0
# 			k ∈ miniFBZ && return exp(-k² / v.α²4[i, i]) / k²
# 			return 0.0
# 		end
# 		(result, err) = hcubature(integrand, a, b; norm = abs, rtol = 1e-8, atol = 1e-10)
# 		head = coff * (√π * v.α[i, i] * erf(qsz / (2 * v.α[i, i])) + result / (4π))
# 		A[i, i] += head
# 	end

# 	return A
# end
function LongRangePart(::Val{:gauss}, lattice::Lattice,
	orblocat::Vector{<:ReducedCoordinates}; norb::Integer = length(orblocat),
	α::Real = 1, αmat::AbstractMatrix{<:Real} = fill(α, (norb, norb)), δ::Real = 1e-8, ϵ::Real = 1)

	Ω = abs(det(parent(lattice)))

	rlattice = reciprocal(lattice)
	rldot = parent(rlattice)
	rldot = Mat3(transpose(rldot) * rldot)

	αmat = _LongRangePartGauss_symmetrize_αmat(αmat)
	φ = [ReciprocalCoulomb(Val(:gauss); ϵ, α = αmat[i, j], Ω) for i in 1:norb, j in 1:norb]
	coff = CoulombScale * 4π / (Ω * ϵ)
	α²4 = 4 .* αmat .* αmat

	Gnorm_max = sqrt(-4 * maximum(αmat)^2 * log(δ)) + norm(parent(rlattice))
	G_frac = gridindex(rlattice, Gnorm_max, Val(3))
	nG = length(G_frac)

	G0idx = findfirst(iszero, G_frac)
	Gn0idx = setdiff(1:nG, G0idx)
	nGminus1 = nG - 1
	G_frac = G_frac[Gn0idx]
	G_car = map(G -> rlattice * G, G_frac)
	G_norm2 = map(G -> G[1] * G[1] + G[2] * G[2] + G[3] * G[3], G_car)

	Δorb = Matrix{ReducedCoordinates{Float64}}(undef, norb, norb)
	Δorb_car = Matrix{CartesianCoordinates{Float64}}(undef, norb, norb)
	G_orb_phase = Array{ComplexF64}(undef, nGminus1, norb, norb)
	Δorb_norm_time2α_divπ = Matrix{Float64}(undef, norb, norb)
	for i in 1:norb, j in 1:i
		dorb = orblocat[i] - orblocat[j]
		dorb_car = lattice * dorb
		Δorb[i, j] = dorb
		Δorb[j, i] = -dorb
		Δorb_car[i, j] = dorb_car
		Δorb_car[j, i] = -dorb_car
		for (Gi, G) in enumerate(G_frac)
			G_orb_phase[Gi, i, j] = cispi(2 * (G ⋅ Δorb[i, j]))
		end
		G_orb_phase[:, j, i] .= conj.(G_orb_phase[:, i, j])
		dorb_norm = norm(dorb_car)
		Δorb_norm_time2α_divπ[i, j] = dorb_norm * 2 * αmat[i, j] / π
		Δorb_norm_time2α_divπ[j, i] = Δorb_norm_time2α_divπ[i, j]
	end

	return LongRangePartGauss(lattice, rlattice, rldot, Ω, norb, orblocat,
		ϵ, αmat, φ, coff, α²4,
		nG, nGminus1, G_frac, G_car, G_norm2,
		Δorb, Δorb_car, G_orb_phase, Δorb_norm_time2α_divπ)
end
function _LongRangePartGauss_symmetrize_αmat(αmat)
	norb = size(αmat, 1)
	αmat_sym = Matrix{Float64}(undef, norb, norb)
	for j in 2:norb, i in 1:(j-1)
		αmat_sym[i, j] = αmat[i, j]
		αmat_sym[j, i] = αmat[i, j]
	end
	for i in 1:norb
		αmat_sym[i, i] = αmat[i, i]
	end
	return αmat_sym
end
