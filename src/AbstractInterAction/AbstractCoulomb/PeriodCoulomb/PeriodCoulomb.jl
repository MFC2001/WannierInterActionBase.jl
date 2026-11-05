export PeriodCoulomb
"""
	PeriodCoulomb <: AbstractCoulomb

Coulomb potential for a periodic system, defined in real space.

It fully accounts for the long-range interactions with all periodic images (mirror systems) 
and removes the divergence via the head-term regularization. 
This implementation can be compared with the result in plane-wave basis.
"""
abstract type PeriodCoulomb <: AbstractCoulomb end

PeriodCoulomb(sym::Symbol, args...; kwargs...) = PeriodCoulomb(Val(sym), args...; kwargs...)

include("./Gauss.jl")
include("./KeldyshGauss.jl")
