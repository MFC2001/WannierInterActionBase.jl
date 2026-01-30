
export MonkhorstPack, RedKgrid, Kline

"""
	abstract type AbstractBrillouinZone end
"""
abstract type AbstractBrillouinZone end

include("./MonkhorstPack.jl")
include("./RedKgrid.jl")
include("./IrredKgrid.jl")
include("./Kline.jl")

"""Bring ``k``-point coordinates into the range [-0.5, 0.5)"""
function normalize_kdirect(x::Real)
	x = x - round(x, RoundNearestTiesUp)
	@assert -1 // 2 ≤ x < 1 // 2
	return x
end
normalize_kdirect(k::AbstractVector{<:Real}) = normalize_kdirect.(k)

"""
Explicitly define the k-points along which to perform BZ sampling.
(Useful for bandstructure calculations)
"""
struct ExplicitKpoints{T} <: AbstractBrillouinZone
	kcoords::Vector{Vec3{T}}
	kweights::Vector{T}
end
function ExplicitKpoints(kcoords::AbstractVector{<:AbstractVector{T}}, kweights::AbstractVector{T}) where {T}
	@assert length(kcoords) == length(kweights)
	ExplicitKpoints{T}(kcoords, kweights)
end
function ExplicitKpoints(kcoords::AbstractVector{<:AbstractVector{T}}) where {T}
	ExplicitKpoints(kcoords, ones(T, length(kcoords)) ./ length(kcoords))
end
function Base.show(io::IO, kgrid::ExplicitKpoints)
	print(io, "ExplicitKpoints with $(length(kgrid.kcoords)) k-points")
end
Base.length(kgrid::ExplicitKpoints) = length(kgrid.kcoords)
