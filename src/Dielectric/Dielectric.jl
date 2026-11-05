export AbstractDielectric

"""
	AbstractDielectric

Abstract supertype of dielectric matrix containers, such as [`DielectricPWBasis`](@ref).
"""
abstract type AbstractDielectric end

include("./DielectricPWBasis.jl")
