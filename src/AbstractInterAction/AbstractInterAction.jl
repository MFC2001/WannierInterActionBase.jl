export AbstractInterAction
"""
	abstract type AbstractInterAction end

Abstract type for different types of interactions. 
"""
abstract type AbstractInterAction end

include("./AbstractCoulomb/AbstractCoulomb.jl")
include("./DirectInterAction/DirectInterAction.jl")
include("./PairInterAction/PairInterAction.jl")
