export PlaneWaveCoulomb
"""
	abstract type PlaneWaveCoulomb <: AbstractCoulomb end

Abstract type for Coulomb interaction potentials under plane-wave basis.

V^{T}(q, G) = \\frac{e^2}{4πΩϵ₀} ∫_{\\mathcal{D}} \\frac{e^{i(q+G)r}}{|r|} dr
"""
abstract type PlaneWaveCoulomb <: AbstractCoulomb end

PlaneWaveCoulomb(sym::Symbol, args...; kwargs...) = PlaneWaveCoulomb(Val(sym), args...; kwargs...)
(v::PlaneWaveCoulomb)(sym::Symbol, args...; kwargs...) = v(Val(sym), args...; kwargs...)

include("./Bulk.jl")
include("./Slab.jl")
