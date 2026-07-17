struct PlaneWaveCoulombBulk <: PlaneWaveCoulomb
	lattice::Lattice{Float64}
	rlattice::ReciprocalLattice{Float64}
	Glist::Vector{ReducedCoordinates{Int}}
	rldot::Mat3{Float64}
	Ω::Float64
	rΩ::Float64
	coff::Float64
end
function (pwc::PlaneWaveCoulombBulk)(q::AbstractVector{<:Real}, Gidx::AbstractVector{<:Integer})
	coulomb = Vector{Float64}(undef, length(Gidx))
	return pwc(coulomb, q, Gidx)
end
function (pwc::PlaneWaveCoulombBulk)(coulomb::AbstractVector{<:Real}, q::AbstractVector{<:Real}, Gidx::AbstractVector{<:Integer})
	q = ReducedCoordinates(q)
	Glist_q = view(pwc.Glist, Gidx)
	Threads.@threads for iG in eachindex(Gidx)
		qG = q + Glist_q[iG]
		qG2 = qG ⋅ (pwc.rldot * qG)
		coulomb[iG] = pwc.coff / qG2
	end
	return coulomb
end
# function (pwc::PlaneWaveCoulombBulk)(::Val{:head}, nq::Integer)
# 	signbit(nq) && error("nq should be positive")
# 	NΩ = nq * pwc.Ω
# 	qsz = (6 * π^2 / NΩ)^(1 // 3)
# 	V_miniFBZ = pwc.rΩ / nq
# 	average = 4π * qsz / V_miniFBZ
# 	return pwc.coff * average
# end
function (pwc::PlaneWaveCoulombBulk)(::Val{:head}, qgrid::Union{AbstractVector, Tuple})
	all(signbit, qgrid) && error("nq should be positive")

	minirlattice = ReciprocalLattice(pwc.rlattice[1:3, 1] ./ qgrid[1], pwc.rlattice[1:3, 2] ./ qgrid[2], pwc.rlattice[1:3, 3] ./ qgrid[3])
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
		q² = q[1] * q[1] + q[2] * q[2] + q[3] * q[3]
		q² < qsz2 && return 0.0
		q ∈ miniFBZ && return 1 / q²
		return 0.0
	end
	(result, err) = hcubature(integrand, min_miniFBZ, max_miniFBZ; norm = abs, rtol = 1e-8, atol = 1e-10)
	result_sphere = 4π * qsz

	nq = prod(qgrid)
	V_miniFBZ = pwc.rΩ / nq

	average = (result + result_sphere) / V_miniFBZ

	return pwc.coff * average
end
function PlaneWaveCoulomb(::Val{:bulk}, lattice::Union{Lattice, AbstractMatrix{<:Real}}, Glist)
	lattice = Lattice(parent(lattice))
	rlattice = reciprocal(lattice)
	rldot = parent(rlattice)
	rldot = Mat3(transpose(rldot) * rldot)
	Ω = abs(det(parent(lattice)))
	rΩ = (2π)^3 / Ω
	Glist = ReducedCoordinates{Int}.(Glist)
	coff = CoulombScale * 4π / Ω
	return PlaneWaveCoulombBulk(lattice, rlattice, Glist, rldot, Ω, rΩ, coff)
end
