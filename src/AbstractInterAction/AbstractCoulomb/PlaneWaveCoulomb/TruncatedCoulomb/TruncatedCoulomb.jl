export TruncatedCoulomb
"""
	abstract type TruncatedCoulomb <: PlaneWaveCoulomb end

Abstract type for truncated Coulomb interaction potentials under plane-wave basis.

V^{T}(q, G) = \\frac{e^2}{4πΩϵ₀} ∫_{\\mathcal{D}} \\frac{e^{i(q+G)r}}{|r|} dr
"""
abstract type TruncatedCoulomb <: PlaneWaveCoulomb end
TruncatedCoulomb(sym::Symbol, args...; kwargs...) = TruncatedCoulomb(Val(sym), args...; kwargs...)

include("./Bare.jl")
include("./Cylinder.jl")
include("./Slab.jl")
include("./Sphere.jl")
