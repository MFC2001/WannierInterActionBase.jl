"""
    HRHermitian{T} <: AbstractHR{T}

A hermitian implementation of `AbstractHR{T}` that stores data of the form `v_{ijR}`.

See the documentation for `AbstractHR` for the meaning of `v_{ijR}`.
"""
struct HRHermitian{T} <: AbstractHR{T}
    n::Int
    R::Vector{ReducedCoordinates{Int}}
    ndegen::Vector{Int}    # If the data is calculated from reciprocal space, then every `R` have a `ndegen`.
    is::Vector{Int}
    js::Vector{Int}
    iRs::Vector{Int}
    vals::Vector{T}
    R2idx::Dict{ReducedCoordinates{Int},Int}
    negR::Vector{Int}
end
@inline _HRHermitian_ispositive(i, j, iR, R) = (i < j || (i == j && _ispositive_vector(R[iR])))
@inline _HRHermitian_isonsite(i, j, iR, R) = (i == j && iszero(R[iR]))
@inline _HRHermitian_isnegative(i, j, iR, R) = (i > j || (i == j && _isnegative_vector(R[iR])))
@inline _HRHermitian_keep(i, j, iR, R) = (i < j || (i == j && !_isnegative_vector(R[iR])))
function Base.show(io::IO, hr::HRHermitian)
    print(io, "Hermitian HR with $(length(hr)) terms and $(numorb(hr)) orbitals.")
end
function _HR_Rindex!(hr::HRHermitian, R::ReducedCoordinates{Int})::Int
    iR = get(hr.R2idx, R, 0)
    if iszero(iR)
        push!(hr.R, R)
        push!(hr.ndegen, 1)
        iR = hr.nR
        hr.R2idx[R] = iR
        push!(hr.R, -R)
        push!(hr.ndegen, 1)
        hr.R2idx[-R] = iR + 1
        push!(hr.negR, iR+1)
        push!(hr.negR, iR)
    end
    return iR
end
function _HR_getindex(hr::HRHermitian{T}, i::Integer, j::Integer, iR::Integer) where {T}
    if i < j
        i_new, j_new, iR_new = i, j, iR
    elseif i == j
        R = hr.R[iR]
        if iszero(R) || _ispositive_vector(R)
            i_new, j_new, iR_new = i, i, iR
        else
            i_new, j_new, iR_new = i, i, hr.negR[iR]
        end
    else
        i_new, j_new, iR_new = j, i, hr.negR[iR]
    end
    idx = findfirst(==((i_new, j_new, iR_new)), indices(hr))
    return isnothing(idx) ? zero(T) : hr.vals[idx]
end
function _HR_setindex!(hr::HRHermitian{T}, v::Number, i::Integer, j::Integer, iR::Integer) where {T}
    if i < j
        i_new, j_new, iR_new, v_new = i, j, iR, v
    elseif i == j
        R = hr.R[iR]
        if iszero(R) || _ispositive_vector(R)
            i_new, j_new, iR_new, v_new = i, i, iR, v
        else
            i_new, j_new, iR_new, v_new = i, i, hr.negR[iR], conj(v)
        end
    else
        i_new, j_new, iR_new, v_new = j, i, hr.negR[iR], conj(v)
    end
    idx = findfirst(==((i_new, j_new, iR_new)), indices(hr))
    isnothing(idx) || return hr.vals[idx] = v_new
    push!(hr.is, i_new)
    push!(hr.js, j_new)
    push!(hr.iRs, iR_new)
    push!(hr.vals, v_new)
    return v_new
end
function _HR_setindex!(hr::HRHermitian{T}, v::Tuple{<:Integer,<:Number}, i::Integer, j::Integer, iR::Integer) where {T}
    if i < j
        i_new, j_new, iR_new, v_new = i, j, iR, v[2]
    elseif i == j
        R = hr.R[iR]
        if iszero(R) || _ispositive_vector(R)
            i_new, j_new, iR_new, v_new = i, i, iR, v[2]
        else
            i_new, j_new, iR_new, v_new = i, i, hr.negR[iR], conj(v[2])
        end
    else
        i_new, j_new, iR_new, v_new = j, i, hr.negR[iR], conj(v[2])
    end
    idx = findfirst(==((i_new, j_new, iR_new)), indices(hr))
    isnothing(idx) || return hr.vals[idx] = v_new
    push!(hr.is, i_new)
    push!(hr.js, j_new)
    push!(hr.iRs, iR_new)
    push!(hr.vals, v_new)
    hr.ndegen[iR] = v[1]
    hr.ndegen[hr.negR[iR]] = v[1]
    return v_new
end
@inline Base.real(hr::HRHermitian{T}) where {T<:Complex} =
    HRHermitian{real(T)}(hr.n, copy(hr.R), copy(hr.ndegen),
        copy(hr.is), copy(hr.js), copy(hr.iRs), real(hr.vals), copy(hr.R2idx), copy(hr.negR))
@inline Base.complex(hr::HRHermitian{T}) where {T<:Real} =
    HRHermitian{complex(T)}(hr.n, copy(hr.R), copy(hr.ndegen),
        copy(hr.is), copy(hr.js), copy(hr.iRs), complex(hr.vals), copy(hr.R2idx), copy(hr.negR))
@inline Base.empty(::Type{HRHermitian{T}}, n::Integer) where {T} =
    HRHermitian{T}(n, ReducedCoordinates{Int}[], Int[], Int[], Int[], Int[], T[], Dict{ReducedCoordinates{Int},Int}(), Int[])
@inline Base.convert(::Type{HRHermitian{T₁}}, hr::HRHermitian{T₂}) where {T₁,T₂} =
    HRHermitian{T₁}(hr.n, copy(hr.R), copy(hr.ndegen),
        copy(hr.is), copy(hr.js), copy(hr.iRs), convert(Vector{T₁}, hr.vals), copy(hr.R2idx), copy(hr.negR))
function Base.filter(f, hr::HRHermitian{T}) where {T}
    mask = collect(f(i, j, iR, v) for (i, j, iR, v) in enumerate(hr))
    is = hr.is[mask]
    js = hr.js[mask]
    iRs = hr.iRs[mask]
    vals = hr.vals[mask]
    return HRHermitian(is, js, iRs, vals, hr.n, hr.R; ndegen=hr.ndegen, inbound=true, dofilter=true, have_negR=true)
end
function Base.sort(hr::HRHermitian{T}) where {T}
    iRs, R, ndegen = _HR_sort_R(hr.iRs, hr.R, hr.ndegen)
    is, js, iRs, vals = _HR_sort_idx(hr.is, hr.js, iRs, hr.vals)
    R2idx = Dict{ReducedCoordinates{Int},Int}(RR=>i for (i, RR) in enumerate(R))
    negR = collect(R2idx[-RR] for RR in R)
    return HRHermitian{T}(hr.n, R, ndegen, is, js, iRs, vals, R2idx, negR)
end
function Base.sort!(hr::HRHermitian{T}) where {T}
    _HR_sort_R!(hr.iRs, hr.R, hr.ndegen)
    _HR_sort_idx!(hr.is, hr.js, hr.iRs, hr.vals)
    for (iR, RR) in enumerate(hr.R)
        hr.R2idx[RR] = iR
    end
    for (iR, RR) in enumerate(hr.R)
        hr.negR[iR] = hr.R2idx[-RR]
    end
    return hr
end
function Base.unique(hr::HRHermitian{T}) where {T}
    unique_path2value = Dict{NTuple{3,Int},T}()
    for (i, j, iR, val) in enumerate(hr)
        path = (i, j, iR)
        if haskey(unique_path2value, path)
            unique_path2value[path] += val
        else
            unique_path2value[path] = val
        end
    end
    n = length(unique_path2value)
    is = Vector{Int}(undef, n)
    js = Vector{Int}(undef, n)
    iRs = Vector{Int}(undef, n)
    vals = Vector{T}(undef, n)
    i = 0
    for (key, value) in unique_path2value
        i += 1
        is[i] = key[1]
        js[i] = key[2]
        iRs[i] = key[3]
        vals[i] = value
    end
    return HRHermitian{T}(hr.n, copy(hr.R), copy(hr.ndegen), is, js, iRs, vals, copy(hr.R2idx), copy(hr.negR))
end
function Base.union(hr₁::HRHermitian{T₁}, hr₂::HRHermitian{T₂}) where {T₁,T₂}
    is = [hr₁.is; hr₂.is .+ numorb(hr₁)]
    js = [hr₁.js; hr₂.js .+ numorb(hr₁)]
    Rs = [hr₁.R[hr₁.iRs]; hr₂.R[hr₂.iRs]]
    vals = [hr₁.vals; hr₂.vals]
    ndegen = [hr₁.ndegen[hr₁.iRs]; hr₂.ndegen[hr₂.iRs]]
    return unique(HRHermitian(is, js, Rs, vals, numorb(hr₁) + numorb(hr₂); ndegen, inbound=true, dofilter=false))
end
function Base.merge(hr₁::HRHermitian{T₁}, hr₂::HRHermitian{T₂}) where {T₁,T₂}
    numorb(hr₁) == numorb(hr₂) || error("Mismatched hr.")
    is = [hr₁.is; hr₂.is]
    js = [hr₁.js; hr₂.js]
    Rs = [hr₁.R[hr₁.iRs]; hr₂.R[hr₂.iRs]]
    vals = [hr₁.vals; hr₂.vals]
    ndegen = [hr₁.ndegen[hr₁.iRs]; hr₂.ndegen[hr₂.iRs]]
    return unique(HRHermitian(is, js, Rs, vals, numorb(hr₁); ndegen, inbound=true, dofilter=false))
end
# 要求输入的数据只有厄密限制下的独立项，否则会重复计数。
function HRHermitian(path::AbstractMatrix{<:Integer}, vals::AbstractVector{T};
    n::Union{Nothing,Integer}=nothing,
    ndegen::Union{Nothing,AbstractVector{<:Integer}}=nothing,
    inbound::Bool=false, dofilter::Bool=true, have_negR::Bool=false,
)::HRHermitian{T} where {T<:Number}

    @assert size(path, 1) == 5 "Wrong size of path."
    if !inbound
        if isnothing(ndegen)
            @assert size(path, 2) == length(vals) "Mismatched length of data."
        else
            @assert all(>(0), ndegen) "ndegen should be positive."
            @assert size(path, 2) == length(vals) == length(ndegen) "Mismatched length of data."
        end
    end

    if isnothing(n)
        length(vals) == 0 && return empty(HRHermitian{T}, 0)
        n = maximum(view(path, 4:5, :))
    else
        length(vals) == 0 && return empty(HRHermitian{T}, n)
        inbound || @assert all(i->1 ≤ i ≤ n, view(path, 4:5, :)) "Invalid is detected: i ∈ [1, $n]."
    end

    is = path[4, :]
    js = path[5, :]
    Rs = map(p->ReducedCoordinates{Int}(p[1], p[2], p[3]), eachcol(path))

    return HRHermitian(is, js, Rs, vals, n; ndegen, inbound=true, dofilter, have_negR)
end
function HRHermitian(
    is::AbstractVector{<:Integer},
    js::AbstractVector{<:Integer},
    Rs::AbstractVector{<:AbstractVector{<:Integer}},
    vals::AbstractVector{T},
    n::Integer;
    ndegen::Union{Nothing,AbstractVector{<:Integer}}=nothing,
    inbound::Bool=false, dofilter::Bool=true, have_negR::Bool=false, docheck::Bool=true,
)::HRHermitian{T} where {T<:Number}

    if !inbound
        @assert all(i->1 ≤ i ≤ n, is) "Invalid is detected: i ∈ [1, $n]."
        @assert all(j->1 ≤ j ≤ n, js) "Invalid js detected: j ∈ [1, $n]."
        if isnothing(ndegen)
            @assert length(is) == length(js) == length(Rs) == length(vals) "Mismatched length of data."
        else
            @assert all(>(0), ndegen) "ndegen should be positive."
            @assert length(is) == length(js) == length(Rs) == length(vals) == length(ndegen) "Mismatched length of data."
        end
    end
    length(vals) == 0 && return empty(HRHermitian{T}, n)

    Rs = collect(ReducedCoordinates{Int}, Rs)
    if dofilter
    end
    R = unique(Rs)
    R = [R; [-RR for RR in R]]
    sort!(R)

    R2idx = Dict{ReducedCoordinates{Int},Int}(RR=>i for (i, RR) in enumerate(R))
    negR = collect(R2idx[-RR] for RR in R)
    iRs = Vector{Int}(undef, length(Rs))
    @inbounds for (i, RR) in enumerate(Rs)
        iRs[i] = R2idx[RR]
    end
    docheck && _HRHermitian_check!(is, js, iRs, vals, R, negR)

    if isnothing(ndegen)
        ndegen_uniqueR = ones(Int, length(R))
    else
        ndegen_uniqueR = Vector{Int}(undef, length(R))
        seen = falses(length(R))
        @inbounds for i in eachindex(iRs)
            iR = iRs[i]
            if seen[iR]
                @assert ndegen_uniqueR[iR] == ndegen[i] "Inconsistent ndegen for the same R vector."
            else
                seen[iR] = true
                ndegen_uniqueR[iR] = ndegen[i]
            end
        end
        for iR in eachindex(R)
            seen[iR] && continue
            ndegen_uniqueR[iR] = ndegen_uniqueR[negR[iR]]
        end
    end
    return HRHermitian{T}(n, R, ndegen_uniqueR, is, js, iRs, vals, R2idx, negR)
end
function HRHermitian(
    is::AbstractVector{<:Integer},
    js::AbstractVector{<:Integer},
    iRs::AbstractVector{<:Integer},
    vals::AbstractVector{T},
    n::Integer,
    R::AbstractVector{<:AbstractVector{<:Integer}};
    ndegen::AbstractVector{<:Integer}=ones(Int, length(R)),
    inbound::Bool=false, dofilter::Bool=true, have_negR::Bool=false, docheck::Bool=true,
)::HRHermitian{T} where {T<:Number}

    if !inbound
        @assert all(i->1 ≤ i ≤ n, is) "Invalid is detected: i ∈ [1, $n]."
        @assert all(j->1 ≤ j ≤ n, js) "Invalid js detected: j ∈ [1, $n]."
        nR = length(R)
        @assert all(iR->1 ≤ iR ≤ nR, iRs) "Invalid ips detected: iR ∈ [1, $nR]."
        @assert length(ndegen) == nR "Mismatched length of ndegen."
        @assert all(>(0), ndegen) "ndegen should be positive."
        @assert length(is) == length(js) == length(iRs) == length(vals) "Mismatched length of data."
    end
    length(vals) == 0 && return empty(HRHermitian{T}, n)

    R = collect(ReducedCoordinates{Int}, R)
    ndegen = collect(Int, ndegen)

    if !have_negR
        iRs, R, ndegen = _HRHermitian_complete_R!(iRs, R, ndegen)
    end
    if dofilter
        # filter R depending on if it exists in iRs, then update iRs and ndegen.
        iRs, R, ndegen = _HRHermitian_filter_R!(iRs, R, ndegen)
    end

    R2idx = Dict{ReducedCoordinates{Int},Int}(RR=>i for (i, RR) in enumerate(R))
    negR = collect(R2idx[-RR] for RR in R)

    docheck && _HRHermitian_check!(is, js, iRs, vals, R, negR)
    return HRHermitian{T}(n, R, ndegen, is, js, iRs, vals, R2idx, negR)
end
# 只存储i<j, i==j && R非负
function _HRHermitian_check!(is, js, iRs, vals, R, negR)
    for ih in eachindex(is)
        i, j, iR = is[ih], js[ih], iRs[ih]
        if _HRHermitian_isnegative(i, j, iR, R)
            is[ih] = j
            js[ih] = i
            iRs[ih] = negR[iR]
            vals[ih] = conj(vals[ih])
        end
    end
end
function _HRHermitian_filter_R!(iRs, R, ndegen)
    R_exist = unique(iRs)
    new_R = unique(view(R, R_exist))
    new_R = unique([new_R; [-RR for RR in new_R]])
    sort!(new_R)
    @inbounds if length(new_R) < length(R)
        new_R_map = Dict{ReducedCoordinates{Int},Int}(RR=>i for (i, RR) in enumerate(new_R))
        index_to_new = [get(new_R_map, RR, 0) for RR in R]
        for (i, iR) in enumerate(iRs)
            iRs[i] = index_to_new[iR]
        end
        R_map = Dict{ReducedCoordinates{Int},Int}(RR=>i for (i, RR) in enumerate(R))
        new_ndegen = [ndegen[R_map[RR]] for RR in new_R]
        return iRs, new_R, new_ndegen
    else
        return iRs, R, ndegen
    end
end
function _HRHermitian_complete_R!(iRs, R, ndegen)
    negR = [-RR for RR in R]
    allR = unique([R; negR])
    sort!(allR)
    allR2idx = Dict{ReducedCoordinates{Int},Int}(RR=>i for (i, RR) in enumerate(allR))
    R_to_allR = [allR2idx[RR] for RR in R]
    negR_to_allR = [allR2idx[RR] for RR in negR]
    allnedegen = zeros(Int, length(allR))
    for (iR, iallR) in enumerate(R_to_allR)
        allnedegen[iallR] = ndegen[iR]
    end
    for (iR, iallR) in enumerate(negR_to_allR)
        if iszero(allnedegen[iallR])
            allnedegen[iallR] = ndegen[iR]
        end
    end
    for (i, iR) in enumerate(iRs)
        iRs[i] = R_to_allR[iR]
    end
    return iRs, allR, allnedegen
end
"""
	spinorize(hr::HRHermitian{T}; mode = conj) -> HRHermitian{T}

Generate a spinful hr from a spinless hr.
Acctually, this function just makes a copy by `mode`, it won't introduce new hopping.
"""
function spinorize(hr::HRHermitian{T}; mode=conj) where {T}
    is = [hr.is; hr.is .+ numorb(hr)]
    js = [hr.js; hr.js .+ numorb(hr)]
    iRs = [hr.iRs; hr.iRs]
    vals = [hr.vals; mode.(hr.vals)]
    return HRHermitian{T}(numorb(hr)*2, copy(hr.R), copy(hr.ndegen), is, js, iRs, vals, copy(hr.R2idx), copy(hr.negR))
end
"""
	translate(hr::HRHermitian{T}, centres::Pair{<:Integer, <:AbstractVector{<:Integer}}...) -> HRHermitian{T}

You only need to input the wannier centres R which isn't in the [0, 0, 0] cell.

See also [`translate!`](@ref).
"""
function translate(hr::HRHermitian{T}, aim_centres::Pair{<:Integer,<:AbstractVector{<:Integer}}...) where {T}
    centres = [ReducedCoordinates{Int}(0, 0, 0) for _ in 1:numorb(hr)]
    for (iw, R) in aim_centres
        centres[iw] = ReducedCoordinates{Int}(R)
    end
    Rs = Vector{ReducedCoordinates{Int}}(undef, length(hr))
    vals = Vector{T}(undef, length(hr))
    for ih in eachindex(Rs)
        i, j, iR = hr.is[ih], hr.js[ih], hr.iRs[ih]
        Rs[ih] = hr.R[iR] + centres[j] - centres[i]
        vals[ih] = hr.vals / hr.ndegen[iR]
    end
    return HRHermitian(copy(hr.is), copy(hr.js), Rs, vals, numorb(hr); inbound=true, dofilter=false)
end
"""
	translate!(hr::HRHermitian, centres::Pair{<:Integer, <:AbstractVector{<:Integer}}...) -> HRHermitian{T}

You only need to input the wannier centres R which isn't in the [0, 0, 0] cell.

See also [`translate`](@ref).
"""
function translate!(hr::HRHermitian{T}, aim_centres::Pair{<:Integer,<:AbstractVector{<:Integer}}...) where {T}
    centres = [ReducedCoordinates{Int}(0, 0, 0) for _ in 1:numorb(hr)]
    for (iw, R) in aim_centres
        centres[iw] = ReducedCoordinates{Int}(R)
    end
    Rs = Vector{ReducedCoordinates{Int}}(undef, length(hr))
    for ih in eachindex(Rs)
        i, j, iR = hr.is[ih], hr.js[ih], hr.iRs[ih]
        Rs[ih] = hr.R[iR] + centres[j] - centres[i]
        hr.vals[ih] /= hr.ndegen[iR]
    end
    empty!(hr.R)
    empty!(hr.R2idx)
    empty!(hr.negR)
    empty!(hr.ndegen)

    allR = unique(Rs)
    allR = unique([allR; [-RR for RR in allR]])

    append!(hr.R, allR)
    sort!(hr.R)
    for (iR, RR) in enumerate(hr.R)
        hr.R2idx[RR] = iR
    end
    append!(hr.negR, collect(hr.R2idx[-RR] for RR in hr.R))
    append!(hr.ndegen, fill(1, hr.nR))
    for (ih, RR) in enumerate(Rs)
        hr.iRs[ih] = hr.R2idx[RR]
    end
    return hr
end
