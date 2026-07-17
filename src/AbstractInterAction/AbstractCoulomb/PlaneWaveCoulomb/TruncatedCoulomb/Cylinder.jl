struct CylinderTruncatedCoulomb <: TruncatedCoulomb
	lattice::Lattice{Float64}
	period::Vec3{Bool}
	Glist::Vector{ReducedCoordinates{Int}}
	R::Float64
	zindex::Int
	xyindex::SVector{2, Int}
	rlatticez::Float64
	Gz::Vector{Int}
	Gxy2::Vector{Float64}
	GRJ1GR::Vector{Float64}
	J0GR::Vector{Float64}
	Ω::Float64
	coff::Float64
end
function (TC::CylinderTruncatedCoulomb)(q::AbstractVector{<:Real}, Gidx::AbstractVector{<:Integer})
	coulomb = Vector{Float64}(undef, length(Gidx))
	return TC(coulomb, q, Gidx)
end
function (TC::CylinderTruncatedCoulomb)(coulomb::AbstractVector{<:Real}, q::AbstractVector{<:Real}, Gidx::AbstractVector{<:Integer})
	qz = q[TC.zindex]
	Gz_q = view(TC.Gz, Gidx)
	Gxy2_q = view(TC.Gxy2, Gidx)
	GRJ1GR_q = view(TC.GRJ1GR, Gidx)
	J0GR_q = view(TC.J0GR, Gidx)
	Threads.@threads for iG in eachindex(Gidx)
		qG_z = TC.rlatticez * (qz + Gz_q[iG])
		qG2 = qG_z * qG_z + Gxy2_q[iG]
		qG_z_R = abs(qG_z) * TC.R
		factor = 1.0 + GRJ1GR_q[iG] * besselk(0, qG_z_R) - qG_z_R * J0GR_q[iG] * besselk(1, qG_z_R)
		coulomb[iG] = TC.coff * factor / qG2
	end
	return coulomb
end
function TruncatedCoulomb(::Val{:cylinder}, lattice::Union{Lattice, AbstractMatrix{<:Real}}, Glist;
	period = [false, false, true], R::Union{Real, Nothing} = nothing)

	check_lattice_period(lattice, period)
	period = Vec3{Bool}(period)
	@assert count(period) == 1 "Wrong period, cylinder truncation is used in 1D system."
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
	rlattice = reciprocal(lattice)

	rlatticez = rlattice[zindex, zindex]
	rlattice_perp = SMatrix{2, 2, Float64, 4}(rlattice[xyindex, xyindex])
	rldot_perp = transpose(rlattice_perp) * rlattice_perp

	if isnothing(R)
		lattice_perp = lattice[xyindex, xyindex]
		h₁ = lattice_perp[1, 1]
		h₂ = lattice_perp[2, 2]
		R = min(h₁, h₂) / 2
	end

	Glist = ReducedCoordinates{Int}.(Glist)
	Gz = map(Glist) do G
		G[zindex]
	end

	nGlist = length(Glist)
	Gxy2 = Vector{Float64}(undef, nGlist)
	GRJ1GR = Vector{Float64}(undef, nGlist)
	J0GR = Vector{Float64}(undef, nGlist)
	Threads.@threads for iG in Base.OneTo(nGlist)
		Gxy = Glist[iG][xyindex]
		Gxy2[iG] = Gxy ⋅ (rldot_perp * Gxy)
		GR = sqrt(Gxy2[iG]) * R
		GRJ1GR[iG] = GR * besselj1(GR)
		J0GR[iG] = besselj0(GR)
	end

	Ω = abs(det(parent(lattice)))
	coff = CoulombScale * 4π / Ω

	return CylinderTruncatedCoulomb(lattice, period, Glist, R,
		zindex, xyindex, rlatticez, Gz, Gxy2, GRJ1GR, J0GR, Ω, coff)
end
