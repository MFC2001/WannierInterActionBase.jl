export AbstractHeadTerm, HeadTerm
"""
	abstract type AbstractHeadTerm <: PlaneWaveCoulomb end

V^{T}(q) = \\frac{e^2}{4πΩϵ₀} ∫_{\\mathcal{D}} \\frac{e^{iqr}}{|r|} dr
\\frac{(2π)^3}{NΩ} \\int_{miniFBZ} V^{T}(q) dk
\\frac{(2π)^2}{NS} \\int_{miniFBZ} V^{T}(q) dk
\\frac{2π}{NL} \\int_{miniFBZ} V^{T}(q) dk
"""
abstract type AbstractHeadTerm <: PlaneWaveCoulomb end
"""
	HeadTerm(sym::Symbol, args...; kwargs...)



Caculate `\\frac{1}{Volume_{miniFBZ}} \\int_{miniFBZ} \\frac{1}{k^2} dk` and so on.
"""
HeadTerm(sym::Symbol, args...; kwargs...) = HeadTerm(Val(sym), args...; kwargs...)


# TODO 调整系数，与TruncatedCoulomb统一


(head::AbstractHeadTerm)() = head.coff * head.value
(head::AbstractHeadTerm)(::Val{:wocoff}) = head.value

include("./Bare.jl")
include("./Box.jl")
include("./Cylinder.jl")
include("./Slab.jl")
include("./Sphere.jl")
include("./Wire.jl")
