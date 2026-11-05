export PairInterAction
"""
	PairInterAction <: AbstractInterAction

Coulomb interaction between pairs defined in reciprocal space,
containing direct-type interaction and pair-type interaction.

\$V_{p_1p_2}(\\mathbf{k})=\\sum_{\\mathbf{R}} V_{p_1\\mathbf{0}p_2\\mathbf{R}} e^{i\\mathbf{k}\\cdot\\mathbf{R}}\$
"""
abstract type PairInterAction <: AbstractInterAction end
function Base.getproperty(int::PairInterAction, name::Symbol)
    if name === :npR
        return length(getfield(int, :pR))
    elseif name === :np
        return length(getfield(int, :pair))
    else
        return getfield(int, name)
    end
end

(int::PairInterAction)(sym::Symbol, args...; kwargs...) = int(Val(sym), args...; kwargs...)

include("./empty.jl")
include("./U.jl")
include("./U_pair.jl")
