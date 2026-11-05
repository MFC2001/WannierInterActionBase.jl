struct PairInterAction_empty <: PairInterAction
    nw::Int
    pR::Vector{ReducedCoordinates{Int}}
    pair::Vector{NTuple{3,Int}}
end
function (int::PairInterAction_empty)(k::AbstractVector{<:Real}, args...; kwargs...)
    return Matrix{ComplexF64}(undef, 0, 0)
end
function (int::PairInterAction_empty)(A::AbstractMatrix{<:Number}, k::AbstractVector{<:Real}, args...; kwargs...)
    size(A) == (0, 0) || error("Buffer size mismatch.")
    return A
end
function PairInterAction(nw::Integer; pair=nothing, atol::Real=0.0)
    pR = ReducedCoordinates{Int}[]
    pU = NTuple{3,Int}[]
    return PairInterAction_empty(nw, pR, pU)
end
