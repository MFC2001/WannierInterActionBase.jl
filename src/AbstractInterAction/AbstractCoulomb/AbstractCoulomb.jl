
"""
	abstract type AbstractCoulomb <: AbstractInterAction end

Abstract type for Coulomb interaction potentials.
"""
abstract type AbstractCoulomb <: AbstractInterAction end

include("./RealCoulomb/RealCoulomb.jl")
include("./ReciprocalCoulomb/ReciprocalCoulomb.jl")
include("./PeriodCoulomb/PeriodCoulomb.jl")
include("./PlaneWaveCoulomb/PlaneWaveCoulomb.jl")
