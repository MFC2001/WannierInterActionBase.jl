export LongRangePart
"""
	abstract type LongRangePart <: AbstractInterAction end

Long-range correction for coulomb interaction, can give the potential matrix in reciprocal space.
It supports 

	(v::LongRangePart)(k::AbstractVector)
	(v::LongRangePart)(k::AbstractVector, nk::Integer)
	(v::LongRangePart)(A::AbstractMatrix, k::AbstractVector)
	(v::LongRangePart)(A::AbstractMatrix, k::AbstractVector, nk::Integer)
	(v::LongRangePart)(::Val{+}, A, k::AbstractVector)
	(v::LongRangePart)(::Val{+}, A, k::AbstractVector, nk::Integer)

The methods involving `nk` can give head term when k is zero in all periodic directions.
"""
abstract type LongRangePart <: AbstractInterAction end

include("./Gauss.jl")
include("./KeldyshGauss.jl")

LongRangePart(sym::Symbol, args...; kwargs...) = LongRangePart(Val(sym), args...; kwargs...)
