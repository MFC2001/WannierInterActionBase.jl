export PeriodCoulomb
"""
	abstract type PeriodCoulomb <: AbstractCoulomb end

The Coulomb potential defined in a periodic system.
"""
abstract type PeriodCoulomb <: AbstractCoulomb end

PeriodCoulomb(sym::Symbol, args...; kwargs...) = PeriodCoulomb(Val(sym), args...; kwargs...)

include("./Gauss.jl")
include("./KeldyshGauss.jl")
