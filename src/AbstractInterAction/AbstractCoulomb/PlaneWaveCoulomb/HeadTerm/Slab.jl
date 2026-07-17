struct SlabHeadTerm <: AbstractHeadTerm
	xyindex::SVector{2, Int}
	zindex::Int
	lattice::SMatrix{2, 2, Float64, 4}
	Rz::Float64
	rlattice::SMatrix{2, 2, Float64, 4}
	qgrid::SVector{2, Int}
	miniFBZ::HalfspaceConvexPolyhedron{Float64, SVector{2, Float64}}
	qsz::Float64
	min_miniFBZ::SVector{2, Float64}
	max_miniFBZ::SVector{2, Float64}
	nq::Int
	S::Float64
	rS::Float64
	Ω::Float64
	value::Float64
	coff::Float64
	α::Float64
end
(head::SlabHeadTerm)(ϵ_∞::Real, q0) = head.coff .* head(Val(:wocoff), ϵ_∞, q0)
function (head::SlabHeadTerm)(::Val{:wocoff}, ϵ_∞::Real, normq0::Real)

	qsz2 = head.qsz * head.qsz
	function integrand(q, γ)
		q² = q[1] * q[1] + q[2] * q[2]
		q² < qsz2 && return 0.0
		q ∈ head.miniFBZ && begin
			q = sqrt(q²)
			return (1 - exp(-q * head.Rz)) / q² / (1 + γ * exp(-head.α * q) * (1 - exp(-q * head.Rz)))
		end
		return 0.0
	end
	integrand_circle(q, γ) = (1 - exp(-q * head.Rz)) / q / (1 + γ * exp(-head.α * q) * (1 - exp(-q * head.Rz)))

	factor_q0 = 1 - exp(-normq0 * head.Rz)

	# v_q0 = (1 - exp(-normq0 * head.Rz)) / normq0^2
	γ = (1 / ϵ_∞ - 1) * exp(head.α * normq0) / factor_q0

	integrand_iω(q) = integrand(q, γ)
	(result, err) = hcubature(integrand_iω, head.min_miniFBZ, head.max_miniFBZ; norm = abs, rtol = 1e-8, atol = 1e-10)

	integrand_circle_iω(q) = integrand_circle(q, γ)
	result_circle, err = quadgk(integrand_circle_iω, 0.0, head.qsz; rtol = 1e-8, atol = 1e-10)

	return (result + 2π * result_circle) * head.nq / head.rS
end
# Construct
function HeadTerm(::Val{:slab}, lattice::Union{Lattice, AbstractMatrix{<:Real}}, qgrid;
	period = [true, true, false], α::Real = 0)

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

	Rz = abs(lattice[zindex, zindex]) / 2
	lattice = lattice[xyindex, xyindex]
	rlattice = rlattice[xyindex, xyindex]

	qgrid = qgrid[xyindex]

	minirlattice = hcat(rlattice[:, 1] ./ qgrid[1], rlattice[:, 2] ./ qgrid[2])
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
			return (1 - exp(-q * Rz)) / q²
		end
		return 0.0
	end
	(result, err) = hcubature(integrand, min_miniFBZ, max_miniFBZ; norm = abs, rtol = 1e-8, atol = 1e-10)

	integrand_circle(q) = (1 - exp(-q * Rz)) / q
	result_circle, err = quadgk(integrand_circle, 0.0, qsz; rtol = 1e-8, atol = 1e-10)

	nq = prod(qgrid)
	𝐚 = lattice[:, 1]
	𝐛 = lattice[:, 2]
	S = abs(𝐚[1] * 𝐛[2] - 𝐚[2] * 𝐛[1])
	𝐚 = rlattice[:, 1]
	𝐛 = rlattice[:, 2]
	rS = abs(𝐚[1] * 𝐛[2] - 𝐚[2] * 𝐛[1])

	value = (result + 2π * result_circle) * nq / rS

	Ω = S * Rz * 2
	coff = 1000 * qₑ / (nq * Ω * ϵ₀)

	return SlabHeadTerm(xyindex, zindex, lattice, Rz, rlattice, qgrid,
		miniFBZ, qsz, min_miniFBZ, max_miniFBZ,
		nq, S, rS, Ω, value, coff, α)
end
