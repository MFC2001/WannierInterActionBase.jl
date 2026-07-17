struct WireHeadTerm <: AbstractHeadTerm
	zindex::Int
	xyindex::SVector{2, Int}
	Rz::Float64
	Rx::Float64
	Ry::Float64
	R::Float64
	rlattice::Float64
	qsz::Float64
	nq::Int
	Ω::Float64
	value::Float64
	coff::Float64
	α::Float64
end
(head::WireHeadTerm)(ϵ_∞::Real, q0) = head.coff .* head(Val(:wocoff), ϵ_∞, q0)
function (head::WireHeadTerm)(::Val{:wocoff}, ϵ_∞::Real, normq0::Real)

	R² = head.R * head.R
	function integrand_vq(q, r)
		r² = r[1] * r[1] + r[2] * r[2]
		r² < R² && return 0.0
		return besselk(0, q * sqrt(r²)) # requier positive q0
	end
	function vq(q)
		(result_q, err) = hcubature(r->integrand_vq(q, r), (0, 0), (head.Rx / 2, head.Ry / 2); norm = abs, rtol = 1e-11, atol = 1e-12)
		result_circle_q = (1 - q * head.R * besselk(1, q * head.R)) * π / (2 * q * q)
		return (result_q + result_circle_q) * 2 / π
	end
	v_q0 = vq(normq0)
	factor_q0 = v_q0 * normq0 * normq0

	function integrand(q, γ)
		vq_value = vq(q)
		ϵq_value = 1 + γ * exp(-head.α * q) * q^2 * vq_value
		return vq_value / ϵq_value
	end

	γ = (1 / ϵ_∞ - 1) * exp(head.α * normq0) / factor_q0

	integrand_iω(q) = integrand(q, γ)
	(result, err) = quadgk(integrand_iω, 0, head.qsz; norm = abs, rtol = 1e-8, atol = 1e-9)

	return result * nq * head.Rz / π
end
# Construct
function HeadTerm(::Val{:wire}, lattice::Union{Lattice, AbstractMatrix{<:Real}}, qgrid;
	period = [false, false, true], α::Real = 0)

	period = Vec3{Bool}(period)
	@assert count(period) == 1 "Wrong period, wire truncation is used in 1D system."
	if period[1]
		zindex = 1
		xyindex = SVector{2, Int}(2, 3)
	elseif period[2]
		zindex = 2
		xyindex = SVector{2, Int}(3, 1)
	else
		zindex = 3
		xyindex = SVector{2, Int}(1, 2)
	end

	lattice = Lattice(parent(lattice))
	Rz = lattice[zindex, zindex]
	Rx = lattice[xyindex[1], xyindex[1]]
	Ry = lattice[xyindex[2], xyindex[2]]
	R = min(Rx, Ry) / 2
	R² = R * R

	rlattice = 2π / Rz

	nq = qgrid[zindex]
	qsz = nq * Rz / π

	function integrand(xyq)
		r² = xyq[1] * xyq[1] + xyq[2] * xyq[2]
		r² < R² && return 0.0
		return besselk(0, xyq[3] * sqrt(r²))
	end
	(result, err) = hcubature(integrand, (0, 0, 0), (Rx/2, Ry/2, qsz); norm = abs, rtol = 1e-8, atol = 1e-10)

	function integrand_cylinder(k)
		kR = k * R
		return 1 - kR * besselk(1, kR) / (k * k)
	end
	(result_cylinder, err) = quadgk(integrand_cylinder, 0, qsz; rtol = 1e-8, atol = 1e-10)

	value = (redult + result_cylinder * π / 2) * 2 * nq * Rz / π^2

	Ω = Rx * Ry * Rz
	coff = 1000 * qₑ / (nq * Ω * ϵ₀)

	return WireHeadTerm(zindex, xyindex, Rz, Rx, Ry, R, rlattice, qsz, nq, Ω, value, coff, α)
end
