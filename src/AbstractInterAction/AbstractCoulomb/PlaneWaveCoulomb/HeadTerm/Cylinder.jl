struct CylinderHeadTerm <: AbstractHeadTerm
	Rz::Float64
	Rx::Float64
	Ry::Float64
	R::Float64
	qsz::Float64
	nq::Int
	Ω::Float64
	value::Float64
	coff::Float64
	α::Float64
end
#TODO How to calculate screened head term?
# Construct
function HeadTerm(::Val{:cylinder}, lattice::Union{Lattice, AbstractMatrix{<:Real}}, qgrid;
	period = [false, false, true], R::Union{Real, Nothing} = nothing, α::Real = 0)

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
	if isnothing(R)
		R = min(Rx, Ry) / 2
	end

	nq = qgrid[zindex]
	qsz = nq * Rz / π

	function integrand(k)
		kR = k * R
		return 1 - kR * besselk(1, kR) / (k * k)
	end
	(result, err) = quadgk(integrand, 0, qsz; rtol = 1e-8, atol = 1e-10)

	value = result / qsz

	Ω = Rx * Ry * Rz
	coff = 1000 * qₑ / (nq * Ω * ϵ₀)

	return CylinderHeadTerm(Rz, Rx, Ry, R, qsz, nq, Ω, value, coff, α)
end
