"""

"""
struct IrredKgrid{UT <: Number} <: AbstractBrillouinZone
	kdirect::Vector{ReducedCoordinates{Rational{Int}}}
	kweight::Vector{Rational{Int}}
	redkdirect::Vector{ReducedCoordinates{Rational{Int}}}
	irmap::Vector{Int}
	symop::Vector{SymOp}
	lattsymop::Vector{LattSymOp{UT}}
	kgrid_size::Vec3{Int}
	kshift::Vec3{Rational{Int}}
end
function Base.show(io::IO, irredkgrid::IrredKgrid)
	print(io, "Irreducibal kgrids with $(length(irredkgrid.kdirect)) irreducible k-points",
		" and $(length(irredkgrid.redkdirect)) reducible k-points.")
end
Base.getindex(irredkgrid::IrredKgrid, index...) = getindex(irredkgrid.kdirect, index...)
Base.length(irredkgrid::IrredKgrid) = length(irredkgrid.kdirect)
Base.iterate(irredkgrid::IrredKgrid, state = 1) = state > length(irredkgrid) ? nothing : (irredkgrid[state], state + 1)
Base.eltype(::IrredKgrid) = ReducedCoordinates{Rational{Int}}
Base.firstindex(::IrredKgrid) = 1
Base.lastindex(irredkgrid::IrredKgrid) = length(irredkgrid)
Base.eachindex(irredkgrid::IrredKgrid) = eachindex(irredkgrid.kdirect)
Base.keys(irredkgrid::IrredKgrid) = keys(irredkgrid.kdirect)
