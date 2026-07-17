struct PlaneWaveCoulombSlab <: PlaneWaveCoulomb
	lattice::Lattice{Float64}
	period::Vec3{Bool}
	Glist::Vector{ReducedCoordinates{Int}}
	α::Float64
	xyindex::SVector{2, Int}
	zindex::Int
	halfRz::Float64
	rlatticexy::SMatrix{2, 2, Float64, 4}
	rldotxy::SMatrix{2, 2, Float64, 4}
	Gxy::Vector{SVector{2, Int}}
	Gz2::Vector{Float64}
	cosGzR::Vector{Bool}
	Ω::Float64
	S::Float64
	rS::Float64
	coff::Float64
end
function (pwc::PlaneWaveCoulombSlab)(q::AbstractVector{<:Real}, Gidx::AbstractVector{<:Integer})
	coulomb = Vector{Float64}(undef, length(Gidx))
	return pwc(coulomb, q, Gidx)
end
function (pwc::PlaneWaveCoulombSlab)(coulomb::AbstractVector{<:Real}, q::AbstractVector{<:Real}, Gidx::AbstractVector{<:Integer})
	qxy = q[pwc.xyindex]
	Gxy_q = view(pwc.Gxy, Gidx)
	Gz2_q = view(pwc.Gz2, Gidx)
	cosGzR_q = view(pwc.cosGzR, Gidx)
	Threads.@threads for iG in eachindex(Gidx)
		qG_xy = qxy + Gxy_q[iG]
		qG_xy2 = qG_xy ⋅ (pwc.rldotxy * qG_xy)
		qG2 = qG_xy2 + Gz2_q[iG]
		factor = cosGzR_q[iG] ? 1.0 - exp(-sqrt(qG_xy2) * pwc.halfRz) : 1.0 + exp(-sqrt(qG_xy2) * pwc.halfRz)
		coulomb[iG] = pwc.coff * factor / qG2
	end
	return coulomb
end
# function (pwc::PlaneWaveCoulombSlab)(::Val{:head}, nq::Integer)
# 	signbit(nq) && error("nq should be positive")
# 	NS = nq * pwc.S
# 	qsz = sqrt(4π / NS)

# 	integrand_circle(q) = (1 - exp(-q * pwc.halfRz)) / q
# 	result_circle, err = quadgk(integrand_circle, 0.0, qsz; rtol = 1e-8, atol = 1e-10)
# 	result_circle *= 2π

# 	S_miniFBZ = pwc.rS / nq
# 	average = result_circle / S_miniFBZ
# 	return pwc.coff * average
# end
function (pwc::PlaneWaveCoulombSlab)(::Val{:head}, qgrid::Union{AbstractVector, Tuple})
	all(signbit, qgrid) && error("nq should be positive")

	qgrid = qgrid[pwc.xyindex]
	minirlattice = hcat(pwc.rlatticexy[:, 1] ./ qgrid[1], pwc.rlatticexy[:, 2] ./ qgrid[2])
	miniFBZ = WignerSeitz(minirlattice)
	qsz2 = minimum(miniFBZ.offsets) / 2
	qsz = sqrt(qsz2)

	min_miniFBZ = Vector{Float64}(undef, 2)
	max_miniFBZ = Vector{Float64}(undef, 2)
	for i in 1:2
		(min_miniFBZ[i], max_miniFBZ[i]) = extrema(k -> k[i], miniFBZ.normals)
	end
	min_miniFBZ = SVector{2, Float64}(min_miniFBZ)
	max_miniFBZ = SVector{2, Float64}(max_miniFBZ)

	function integrand(q)
		q² = q[1]^2 + q[2]^2
		q² < qsz2 && return 0.0
		q ∈ miniFBZ && begin
			q = sqrt(q²)
			return (1 - exp(-q * pwc.halfRz)) / q²
		end
		return 0.0
	end
	(result, err) = hcubature(integrand, min_miniFBZ, max_miniFBZ; norm = abs, rtol = 1e-8, atol = 1e-10)

	integrand_circle(q) = (1 - exp(-q * pwc.halfRz)) / q
	result_circle, err = quadgk(integrand_circle, 0.0, qsz; rtol = 1e-8, atol = 1e-10)
	result_circle *= 2π

	nq = prod(qgrid)
	S_miniFBZ = pwc.rS / nq

	average = (result + result_circle) / S_miniFBZ
	return pwc.coff * average
end
function (pwc::PlaneWaveCoulombSlab)(::Val{:head}, qgrid::Union{AbstractVector, Tuple}, ϵ_∞::Real, normq0::Real)
	all(signbit, qgrid) && error("nq should be positive")

	qgrid = qgrid[pwc.xyindex]
	minirlattice = hcat(pwc.rlatticexy[:, 1] ./ qgrid[1], pwc.rlatticexy[:, 2] ./ qgrid[2])
	miniFBZ = WignerSeitz(minirlattice)
	qsz2 = minimum(miniFBZ.offsets) / 2
	qsz = sqrt(qsz2)

	min_miniFBZ = Vector{Float64}(undef, 2)
	max_miniFBZ = Vector{Float64}(undef, 2)
	for i in 1:2
		(min_miniFBZ[i], max_miniFBZ[i]) = extrema(k -> k[i], miniFBZ.normals)
	end
	min_miniFBZ = SVector{2, Float64}(min_miniFBZ)
	max_miniFBZ = SVector{2, Float64}(max_miniFBZ)

	# γ
	factor_q0 = 1 - exp(-normq0 * pwc.halfRz)
	# v_q0 = (1 - exp(-normq0 * pwc.halfRz)) / normq0^2
	γ = (1 / ϵ_∞ - 1) * exp(pwc.α * normq0) / factor_q0

	function integrand(q)
		q² = q[1] * q[1] + q[2] * q[2]
		q² < qsz2 && return 0.0
		q ∈ miniFBZ && begin
			q = sqrt(q²)
			return (1 - exp(-q * pwc.halfRz)) / q² / (1 + γ * exp(-pwc.α * q) * (1 - exp(-q * pwc.halfRz)))
		end
		return 0.0
	end
	(result, err) = hcubature(integrand, min_miniFBZ, max_miniFBZ; norm = abs, rtol = 1e-8, atol = 1e-10)

	integrand_circle(q) = (1 - exp(-q * pwc.halfRz)) / q / (1 + γ * exp(-pwc.α * q) * (1 - exp(-q * pwc.halfRz)))
	(result_circle, err) = quadgk(integrand_circle, 0.0, qsz; rtol = 1e-8, atol = 1e-10)
	result_circle *= 2π

	nq = prod(qgrid)
	S_miniFBZ = pwc.rS / nq

	average = (result + result_circle) / S_miniFBZ
	return pwc.coff * average
end
function (pwc::PlaneWaveCoulombSlab)(::Val{:head}, qgrid::Union{AbstractVector, Tuple}, ϵ_∞::AbstractVector{<:Real}, normq0::Real)
	all(signbit, qgrid) && error("nq should be positive")

	qgrid = qgrid[pwc.xyindex]
	minirlattice = hcat(pwc.rlatticexy[:, 1] ./ qgrid[1], pwc.rlatticexy[:, 2] ./ qgrid[2])
	miniFBZ = WignerSeitz(minirlattice)
	qsz2 = minimum(miniFBZ.offsets) / 2
	qsz = sqrt(qsz2)

	min_miniFBZ = Vector{Float64}(undef, 2)
	max_miniFBZ = Vector{Float64}(undef, 2)
	for i in 1:2
		(min_miniFBZ[i], max_miniFBZ[i]) = extrema(k -> k[i], miniFBZ.normals)
	end
	min_miniFBZ = SVector{2, Float64}(min_miniFBZ)
	max_miniFBZ = SVector{2, Float64}(max_miniFBZ)

	function integrand(q, γ)
		q² = q[1] * q[1] + q[2] * q[2]
		q² < qsz2 && return 0.0
		q ∈ miniFBZ && begin
			q = sqrt(q²)
			return (1 - exp(-q * pwc.halfRz)) / q² / (1 + γ * exp(-pwc.α * q) * (1 - exp(-q * pwc.halfRz)))
		end
		return 0.0
	end
	integrand_circle(q, γ) = (1 - exp(-q * pwc.halfRz)) / q / (1 + γ * exp(-pwc.α * q) * (1 - exp(-q * pwc.halfRz)))

	nω = length(ϵ_∞)
	average = Vector{Float64}(undef, nω)
	for iω in Base.OneTo(nω)
		# γ
		factor_q0 = 1 - exp(-normq0 * pwc.halfRz)
		# v_q0 = (1 - exp(-normq0 * pwc.halfRz)) / normq0^2
		γ = (1 / ϵ_∞[iω] - 1) * exp(pwc.α * normq0) / factor_q0

		integrand_iω(q) = integrand(q, γ)
		(result, err) = hcubature(integrand_iω, min_miniFBZ, max_miniFBZ; norm = abs, rtol = 1e-8, atol = 1e-10)
		integrand_circle_iω(q) = integrand_circle(q, γ)
		(result_circle, err) = quadgk(integrand_circle_iω, 0.0, qsz; rtol = 1e-8, atol = 1e-10)
		result_circle *= 2π

		average[iω] = result + result_circle
	end

	nq = prod(qgrid)
	S_miniFBZ = pwc.rS / nq
	average ./= S_miniFBZ

	return pwc.coff .* average
end
function PlaneWaveCoulomb(::Val{:slab}, lattice::Union{Lattice, AbstractMatrix{<:Real}}, Glist;
	period = [true, true, false], α::Real = 0)
	check_lattice_period(lattice, period)
	period = Vec3{Bool}(period)
	@assert count(period) == 2 "Wrong period, slab truncation is used in 2D system."
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

	lattice = Lattice(parent(lattice))
	rlattice = reciprocal(lattice)

	latticexy = SMatrix{2, 2, Float64, 4}(lattice[xyindex, xyindex])
	halfRz = abs(lattice[zindex, zindex]) / 2
	rlatticexy = SMatrix{2, 2, Float64, 4}(rlattice[xyindex, xyindex])
	rldotxy = transpose(rlatticexy) * rlatticexy
	rldotz = rlattice[zindex, zindex] * rlattice[zindex, zindex]

	Glist = ReducedCoordinates{Int}.(Glist)
	Gxy = map(Glist) do G
		G[xyindex]
	end
	Gz2 = map(Glist) do G
		Gz = G[zindex]
		Gz * rldotz * Gz
	end
	cosGzR = map(Glist) do G
		iseven(G[zindex]) ? true : false
	end

	𝐚 = latticexy[:, 1]
	𝐛 = latticexy[:, 2]
	S = abs(𝐚[1] * 𝐛[2] - 𝐚[2] * 𝐛[1])
	𝐚 = rlatticexy[:, 1]
	𝐛 = rlatticexy[:, 2]
	rS = abs(𝐚[1] * 𝐛[2] - 𝐚[2] * 𝐛[1])

	Ω = S * halfRz * 2
	coff = CoulombScale * 4π / Ω

	return PlaneWaveCoulombSlab(lattice, period, Glist, α,
		xyindex, zindex, halfRz, rlatticexy, rldotxy, Gxy, Gz2, cosGzR, Ω, S, rS, coff)
end
