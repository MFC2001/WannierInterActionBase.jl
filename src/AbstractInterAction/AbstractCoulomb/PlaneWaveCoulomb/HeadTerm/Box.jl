struct BoxHeadTerm <: AbstractHeadTerm
	lattice::Lattice{Float64}
	R::Float64
	Ω::Float64
	value::Float64
	coff::Float64
end
#TODO How to calculate screened head term?
# Construct
function HeadTerm(::Val{:box}, lattice::Union{Lattice, AbstractMatrix{<:Real}})
	lattice = Lattice(parent(lattice))
	Ω = abs(det(parent(lattice)))
	h₁ = lattice[1, 1]
	h₂ = lattice[2, 2]
	h₃ = lattice[3, 3]
	R = min(h₁, h₂, h₃) / 2

	min_unitcell = Vec3(-h₁/2, -h₂/2, -h₃/2)
	max_unitcell = Vec3(h₁/2, h₂/2, h₃/2)

	R² = R * R
	function integrand(r)
		r² = r[1]^2 + r[2]^2 + r[3]^2
		r² < R² && return 0.0
		return 1 / sqrt(r²)
	end

	(result, err) = hcubature(integrand, min_unitcell, max_unitcell; norm = abs, rtol = 1e-8, atol = 1e-10)

	value = (result + 2π * R²) / 4π
	coff = 1000 * qₑ / (Ω * ϵ₀)

	return BoxHeadTerm(lattice, R, Ω, value, coff)
end
