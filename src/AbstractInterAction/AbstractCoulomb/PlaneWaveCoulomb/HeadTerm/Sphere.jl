struct SphereHeadTerm <: AbstractHeadTerm
	lattice::Lattice{Float64}
	R::Float64
	Ω::Float64
	value::Float64
	coff::Float64
end
#TODO How to calculate screened head term?
# Construct
function HeadTerm(::Val{:sphere}, lattice::Union{Lattice, AbstractMatrix{<:Real}}; R::Union{Real, Nothing} = nothing)
	lattice = Lattice(parent(lattice))
	Ω = abs(det(parent(lattice)))
	if isnothing(R)
		h₁ = lattice[1, 1]
		h₂ = lattice[2, 2]
		h₃ = lattice[3, 3]
		R = min(h₁, h₂, h₃) / 2
	end

	value = R * R / 2
	coff = 1000 * qₑ / (Ω * ϵ₀)

	return SphereHeadTerm(lattice, R, Ω, value, coff)
end
