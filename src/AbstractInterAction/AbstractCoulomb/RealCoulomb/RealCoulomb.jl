export RealCoulomb
"""
	RealCoulomb <: AbstractCoulomb

The Coulomb potential defined in real space.
"""
abstract type RealCoulomb <: AbstractCoulomb end

RealCoulomb(sym::Symbol, args...; kwargs...) = RealCoulomb(Val(sym), args...; kwargs...)

include("./Erfc.jl")
include("./Gauss.jl")
include("./Keldysh.jl")
include("./KeldyshGauss.jl")
include("./Ohno.jl")
include("./Point.jl")
include("./SoftCore.jl")
