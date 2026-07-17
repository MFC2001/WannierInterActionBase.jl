struct SphereTruncatedCoulomb <: TruncatedCoulomb
	lattice::Lattice{Float64}
	Glist::Vector{ReducedCoordinates{Int}}
	R::Float64
	coulomb::Vector{Float64}
end
(TC::SphereTruncatedCoulomb)() = TC.coulomb
function TruncatedCoulomb(::Val{:sphere}, lattice::Union{Lattice, AbstractMatrix{<:Real}}, Glist; R::Union{Real, Nothing} = nothing)
	check_lattice_period(lattice, [false, false, false])
	if isnothing(R)
		h₁ = lattice[1, 1]
		h₂ = lattice[2, 2]
		h₃ = lattice[3, 3]
		R = min(h₁, h₂, h₃) / 2
	end

	lattice = Lattice(parent(lattice))
	rlattice = reciprocal(lattice)

	Ω = abs(det(parent(lattice)))
	coff = CoulombScale * 4π / Ω

	rlx = rlattice[1, 1]
	rly = rlattice[2, 2]
	rlz = rlattice[3, 3]
	rlx2 = rlx * rlx
	rly2 = rly * rly
	rlz2 = rlz * rlz

	Glist = ReducedCoordinates{Int}.(Glist)
	coulomb = map(Glist) do G
		G2 = G[1] * rlx2 * G[1] + G[2] * rly2 * G[2] + G[3] * rlz2 * G[3]
		normG = sqrt(G2)
		coff * (1.0 - cos(normG * R)) / G2
	end
	return SphereTruncatedCoulomb(lattice, Glist, R, coulomb)
end
