export RealCoulomb
"""
	abstract type RealCoulomb <: AbstractCoulomb end

The Coulomb potential defined in real space.
"""
abstract type RealCoulomb <: AbstractCoulomb end

RealCoulomb(sym::Symbol, args...; kwargs...) = RealCoulomb(Val(sym), args...; kwargs...)

include("./Point.jl")
include("./Gauss.jl")
include("./Keldysh.jl")
include("./KeldyshGauss.jl")
