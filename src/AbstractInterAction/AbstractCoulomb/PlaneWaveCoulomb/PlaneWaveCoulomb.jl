export PlaneWaveCoulomb
"""
	PlaneWaveCoulomb <: AbstractCoulomb

The Coulomb potential defined in reciprocal space, 
describing the Coulomb interaction in plane-wave basis.

V^{T}(q, G) = \\frac{e^2}{4πΩϵ₀} ∫_{\\mathcal{D}} \\frac{e^{i(q+G)r}}{|r|} dr

Note that the 1D and 0D branches have not been tested yet.
"""
abstract type PlaneWaveCoulomb <: AbstractCoulomb end

PlaneWaveCoulomb(sym::Symbol, args...; kwargs...) = PlaneWaveCoulomb(Val(sym), args...; kwargs...)
(v::PlaneWaveCoulomb)(sym::Symbol, args...; kwargs...) = v(Val(sym), args...; kwargs...)

include("./Box.jl")
include("./Bulk.jl")
include("./Cylinder.jl")
include("./Slab.jl")
include("./Sphere.jl")
include("./Wire.jl")
