export IdealCell, IdealHR

"""
    IdealCell(sym::Symbol, args...; kwargs...)

See details in `/IdealLattice`.
"""
IdealCell(sym::Symbol, args...; kwargs...) = IdealCell(Val(sym), args...; kwargs...)
"""
    IdealHR(sym::Symbol, args...; kwargs...)

See details in `/IdealLattice`.
"""
IdealHR(sym::Symbol, args...; kwargs...) = IdealHR(Val(sym), args...; kwargs...)

include("./gammagraphyne.jl")
include("./hexagonal.jl")
include("./kagome_hexagonal.jl")
include("./kagome.jl")
include("./star.jl")
