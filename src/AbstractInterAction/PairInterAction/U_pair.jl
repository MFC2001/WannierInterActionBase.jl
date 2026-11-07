struct PairInterAction_U_pair{directT<:DirectInterAction,pairT<:AbstractReciprocalHoppings} <: PairInterAction
    nw::Int
    pR::Vector{ReducedCoordinates{Int}}
    pair::Vector{NTuple{3,Int}}
    U_index::Vector{Int}
    Uint::directT
    pairint::pairT
end
function (int::PairInterAction_U_pair)(k::AbstractVector{<:Real}, args...; kwargs...)
    A = int.pairint(k)
    Umat = int.Uint(k, args...; kwargs...)
    A[int.U_index, int.U_index] .= Umat
    return A
end
function (int::PairInterAction_U_pair)(A::AbstractMatrix{<:Number}, k::AbstractVector{<:Real}, args...; kwargs...)
    int.pairint(A, k)
    Umat = int.Uint(k, args...; kwargs...)
    A[int.U_index, int.U_index] .= Umat
    return A
end
"""
    PairInterAction(Uint::DirectInterAction, pairint::PairIntHR; pair=nothing, atol::Real=0.0)

Construct an object containing all interactions.

- `pair`: a vector of `(i, j, R)`, where `i` and `j` are Wannier indices and `R` is a 3-vector of integers. If given, only interactions for these pairs are retained. See also [`filter_pair`](@ref).
- `atol`: entries in `pairint` with absolute value below `atol` are discarded. See also [`prune`](@ref).
"""
function PairInterAction(Uint::DirectInterAction, pairint::PairIntHRFull{T}; pair=nothing, atol::Real=0.0) where {T}
    nw = numorb(Uint)
    @assert numorb(pairint) == nw "Mismatched direct-type and pair-type interaction."
    if !isnothing(pair)
        pairint = filter_pair(pairint, pair)
    end
    if atol > 0
        pairint = prune(pairint, atol)
    end

    pairint_pair = map(p -> (p[1], p[2], pairint.pR[p[3]]), pairint.pair)
    pU = [(i, i, ReducedCoordinates{Int}(0, 0, 0)) for i in 1:nw]
    allpair = [pU; unique(setdiff(pairint_pair, pU))]
    nallpair = length(allpair)
    allpair_idx = Dict{Tuple{Int,Int,ReducedCoordinates{Int}},Int}(p=>i for (i, p) in enumerate(allpair))

    U_index = map(p->allpair_idx[p], pU)
    pair_index = map(p->allpair_idx[p], pairint_pair)

    hops = [Vector{_Hopping_ij{Int,T}}(undef, 0) for _ in CartesianIndices((nallpair, nallpair))]
    nR = pairint.nR
    for jp in Base.OneTo(nallpair), ip in Base.OneTo(nallpair)
        sizehint!(hops[ip, jp], nR)
    end
    for (i, j, iR, val) in enumerate(pairint)
        i_new, j_new = pair_index[i], pair_index[j]
        push!(hops[i_new, j_new], _Hopping_ij{Int,T}(pairint.R[iR], val/pairint.ndegen[iR]))
    end
    Nhop = length.(hops)
    pairint = BaseReciprocalHoppings{T}(nallpair, Nhop, hops)
    pairint = ReciprocalHoppings(pairint)

    allpR = unique([p[3] for p in allpair])
    allpR_idx = Dict{ReducedCoordinates{Int},Int}(pR=>i for (i, pR) in enumerate(allpR))
    allpair = map(p->(p[1], p[2], allpR_idx[p[3]]), allpair)

    return PairInterAction_U_pair{typeof(Uint),typeof(pairint)}(nw, allpR, allpair, U_index, Uint, pairint)
end
function PairInterAction(Uint::DirectInterAction, pairint::PairIntHRHermitian{T}; pair=nothing, atol::Real=0.0) where {T}
    nw = numorb(Uint)
    @assert numorb(pairint) == nw "Mismatched direct-type and pair-type interaction."
    if !isnothing(pair)
        pairint = filter_pair(pairint, pair)
    end
    if atol > 0
        pairint = prune(pairint, atol)
    end

    pairint_pair = map(p -> (p[1], p[2], pairint.pR[p[3]]), pairint.pair)
    pU = [(i, i, ReducedCoordinates{Int}(0, 0, 0)) for i in 1:nw]
    allpair = [pU; unique(setdiff(pairint_pair, pU))]
    nallpair = length(allpair)
    allpair_idx = Dict{Tuple{Int,Int,ReducedCoordinates{Int}},Int}(p=>i for (i, p) in enumerate(allpair))

    U_index = map(p->allpair_idx[p], pU)
    pair_index = map(p->allpair_idx[p], pairint_pair)

    hops = [Vector{_Hopping_ij{Int,T}}(undef, 0) for _ in CartesianIndices((nallpair, nallpair))]
    nR = pairint.nR
    for jp in Base.OneTo(nallpair), ip in 1:jp
        sizehint!(hops[ip, jp], nR)
    end
    notzero = map(x->!iszero(x), pairint.R)
    for (i, j, iR, val) in enumerate(pairint)
        i_new, j_new = pair_index[i], pair_index[j]
        # new index for allpair, i < j ≠ i_new < j_new
        if i_new > j_new
            iR = pairint.negR[iR]
            push!(hops[j_new, i_new], _Hopping_ij{Int,T}(pairint.R[iR], conj(val/pairint.ndegen[iR])))
        else
            push!(hops[i_new, j_new], _Hopping_ij{Int,T}(pairint.R[iR], val/pairint.ndegen[iR]))
            if i_new == j_new && notzero[iR]
                iR = pairint.negR[iR]
                push!(hops[i_new, i_new], _Hopping_ij{Int,T}(pairint.R[iR], conj(val/pairint.ndegen[iR])))
            end
        end
    end
    Nhop = length.(hops)
    pairint = BaseReciprocalHoppings{T}(nallpair, Nhop, hops)
    pairint = HermitianReciprocalHoppings(pairint, :U)

    allpR = unique([p[3] for p in allpair])
    allpR_idx = Dict{ReducedCoordinates{Int},Int}(pR=>i for (i, pR) in enumerate(allpR))
    allpair = map(p->(p[1], p[2], allpR_idx[p[3]]), allpair)

    return PairInterAction_U_pair{typeof(Uint),typeof(pairint)}(nw, allpR, allpair, U_index, Uint, pairint)
end
