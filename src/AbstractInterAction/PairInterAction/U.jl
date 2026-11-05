struct PairInterAction_U{directT<:DirectInterAction} <: PairInterAction
    nw::Int
    pR::Vector{ReducedCoordinates{Int}}
    pair::Vector{NTuple{3,Int}}
    Uint::directT
end
function (int::PairInterAction_U)(k::AbstractVector{<:Real}, args...; kwargs...)
    return int.Uint(k, args...; kwargs...)
end
function (int::PairInterAction_U)(A::AbstractMatrix{<:Number}, k::AbstractVector{<:Real}, args...; kwargs...)
    return int.Uint(A, k, args...; kwargs...)
end
function PairInterAction(Uint::DirectInterAction; pair=nothing, atol::Real=0.0)
    nw = numorb(Uint)
    pR = [ReducedCoordinates{Int}(0, 0, 0)]
    pU = [(i, i, 1) for i in 1:nw]
    return PairInterAction_U{typeof(Uint)}(nw, pR, pU, Uint)
end
