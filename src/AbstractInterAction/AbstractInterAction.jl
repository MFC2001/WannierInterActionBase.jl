export AbstractInterAction
"""
	abstract type AbstractInterAction end

Abstract type for interaction potentials. 
"""
abstract type AbstractInterAction end

include("./AbstractCoulomb/AbstractCoulomb.jl")
include("./DirectInterAction/DirectInterAction.jl")
include("./MirrorCorrection/MirrorCorrection.jl")
include("./LRcheck.jl")
