struct BareHeadTerm <: AbstractHeadTerm
	lattice::Lattice{Float64}
	rlattice::ReciprocalLattice{Float64}
	qgrid::Vec3{Int}
	miniFBZ::HalfspaceConvexPolyhedron{Float64, CartesianCoordinates{Float64}}
	qsz::Float64
	min_miniFBZ::Vec3{Float64}
	max_miniFBZ::Vec3{Float64}
	nq::Int
	Ω::Float64
	rΩ::Float64
	value::Float64
	coff::Float64
end
# (head::BareHeadTerm)(ϵ_∞, q0) = head.coff .* head(Val(:wocoff), ϵ_∞)
# (head::BareHeadTerm)(::Val{:wocoff}, ϵ_∞, q0) = head(Val(:wocoff), ϵ_∞)
# (head::BareHeadTerm)(ϵ_∞) = head.coff .* head(Val(:wocoff), ϵ_∞)
# function (head::BareHeadTerm)(::Val{:wocoff}, ϵ_∞)
# 	nω = length(ϵ_∞)
# 	head_nω1 = Vector{ComplexF64}(undef, nω + 1)
# 	for iω in 1:nω
# 		head_nω1[iω] = head.value * ϵ_∞[iω]
# 	end
# 	head_nω1[end] = head.value
# 	return head_nω1
# end
# Construct
function HeadTerm(::Val{:bare}, lattice::Union{Lattice, AbstractMatrix{<:Real}}, qgrid)
	lattice = Lattice(parent(lattice))
	qgrid = Vec3{Int}(qgrid[1], qgrid[2], qgrid[3])
	rlattice = reciprocal(lattice)

	minirlattice = ReciprocalLattice(rlattice[1:3, 1] ./ qgrid[1], rlattice[1:3, 2] ./ qgrid[2], rlattice[1:3, 3] ./ qgrid[3])
	miniFBZ = WignerSeitz(minirlattice)

	qsz2 = minimum(miniFBZ.offsets) / 2
	qsz = sqrt(qsz2)

	min_miniFBZ = Vector{Float64}(undef, 3)
	max_miniFBZ = Vector{Float64}(undef, 3)
	for i in 1:3
		(min_miniFBZ[i], max_miniFBZ[i]) = extrema(k -> k[i], miniFBZ.normals)
	end
	min_miniFBZ = Vec3(min_miniFBZ)
	max_miniFBZ = Vec3(max_miniFBZ)

	function integrand(q)
		q² = q[1]^2 + q[2]^2 + q[3]^2
		q² < qsz2 && return 0.0
		q ∈ miniFBZ && return 1 / q²
		return 0.0
	end

	(result, err) = hcubature(integrand, min_miniFBZ, max_miniFBZ; norm = abs, rtol = 1e-8, atol = 1e-10)

	nq = prod(qgrid)

	Ω = abs(det(parent(lattice)))
	rΩ = abs(det(parent(rlattice)))

	value = (result + 4π * qsz) * nq / rΩ

	coff = CoulombScale * 4π / Ω

	return BareHeadTerm(lattice, rlattice, qgrid,
		miniFBZ, qsz, min_miniFBZ, max_miniFBZ,
		nq, Ω, rΩ, value, coff)
end
