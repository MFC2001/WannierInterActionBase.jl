"""
    HRFull{T} <: AbstractHR{T}

A simple implementation of `AbstractHR{T}` that stores data of the form `v_{ijR}`.

See the documentation for `AbstractHR` for the meaning of `v_{ijR}`.
"""
struct HRFull{T} <: AbstractHR{T}
    n::Int
    R::Vector{ReducedCoordinates{Int}}
    ndegen::Vector{Int}    # If the data is calculated from reciprocal space, then every `R` have a `ndegen`.
    is::Vector{Int}
    js::Vector{Int}
    iRs::Vector{Int}
    vals::Vector{T}
    R2idx::Dict{ReducedCoordinates{Int},Int}
end
function Base.show(io::IO, hr::HRFull)
    print(io, "HR with $(length(hr)) terms and $(numorb(hr)) orbitals.")
end
function _HR_Rindex!(hr::HRFull, R::ReducedCoordinates{Int})::Int
    iR = get(hr.R2idx, R, 0)
    if iszero(iR)
        push!(hr.R, R)
        push!(hr.ndegen, 1)
        iR = hr.nR
        hr.R2idx[R] = iR
    end
    return iR
end
function _HR_getindex(hr::HRFull{T}, i::Integer, j::Integer, iR::Integer) where {T}
    idx = findfirst(==((i, j, iR)), indices(hr))
    return isnothing(idx) ? zero(T) : hr.vals[idx]
end
function _HR_setindex!(hr::HRFull{T}, v::Number, i::Integer, j::Integer, iR::Integer) where {T}
    idx = findfirst(==((i, j, iR)), indices(hr))
    isnothing(idx) || return hr.vals[idx] = v
    push!(hr.is, i)
    push!(hr.js, j)
    push!(hr.iRs, iR)
    push!(hr.vals, v)
    return v
end
function _HR_setindex!(hr::HRFull{T}, v::Tuple{<:Integer,<:Number}, i::Integer, j::Integer, iR::Integer) where {T}
    idx = findfirst(==((i, j, iR)), indices(hr))
    isnothing(idx) || return hr.vals[idx] = v[2]
    push!(hr.is, i)
    push!(hr.js, j)
    push!(hr.iRs, iR)
    push!(hr.vals, v[2])
    hr.ndegen[iR] = v[1]
    return v[2]
end
@inline Base.real(hr::HRFull{T}) where {T<:Complex} =
    HRFull{real(T)}(hr.n, copy(hr.R), copy(hr.ndegen),
        copy(hr.is), copy(hr.js), copy(hr.iRs), real(hr.vals), copy(hr.R2idx))
@inline Base.complex(hr::HRFull{T}) where {T<:Real} =
    HRFull{complex(T)}(hr.n, copy(hr.R), copy(hr.ndegen),
        copy(hr.is), copy(hr.js), copy(hr.iRs), complex(hr.vals), copy(hr.R2idx))
@inline Base.empty(::Type{HRFull{T}}, n::Integer) where {T} =
    HRFull{T}(n, ReducedCoordinates{Int}[], Int[], Int[], Int[], Int[], T[], Dict{ReducedCoordinates{Int},Int}())
@inline Base.convert(::Type{HRFull{T₁}}, hr::HRFull{T₂}) where {T₁,T₂} =
    HRFull{T₁}(hr.n, copy(hr.R), copy(hr.ndegen),
        copy(hr.is), copy(hr.js), copy(hr.iRs), convert(Vector{T₁}, hr.vals), copy(hr.R2idx))
function Base.filter(f, hr::HRFull{T}) where {T}
    mask = collect(f(i, j, iR, v) for (i, j, iR, v) in enumerate(hr))
    is = hr.is[mask]
    js = hr.js[mask]
    iRs = hr.iRs[mask]
    vals = hr.vals[mask]
    return HRFull(is, js, iRs, vals, hr.n, hr.R; ndegen=hr.ndegen, inbound=true, dofilter=true)
end
function Base.sort(hr::HRFull{T}) where {T}
    iRs, R, ndegen = _HR_sort_R(hr.iRs, hr.R, hr.ndegen)
    is, js, iRs, vals = _HR_sort_idx(hr.is, hr.js, iRs, hr.vals)
    R2idx = Dict{ReducedCoordinates{Int},Int}(RR=>i for (i, RR) in enumerate(R))
    return HRFull{T}(hr.n, R, ndegen, is, js, iRs, vals, R2idx)
end
function Base.sort!(hr::HRFull{T}) where {T}
    _HR_sort_R!(hr.iRs, hr.R, hr.ndegen)
    _HR_sort_idx!(hr.is, hr.js, hr.iRs, hr.vals)
    for (iR, RR) in enumerate(hr.R)
        hr.R2idx[RR] = iR
    end
    return hr
end
function Base.unique(hr::HRFull{T}) where {T}
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
    return HRFull{T}(hr.n, copy(hr.R), copy(hr.ndegen), is, js, iRs, vals, copy(hr.R2idx))
end
function Base.union(hr₁::HRFull{T₁}, hr₂::HRFull{T₂}) where {T₁,T₂}
    is = [hr₁.is; hr₂.is .+ numorb(hr₁)]
    js = [hr₁.js; hr₂.js .+ numorb(hr₁)]
    Rs = [hr₁.R[hr₁.iRs]; hr₂.R[hr₂.iRs]]
    vals = [hr₁.vals; hr₂.vals]
    ndegen = [hr₁.ndegen[hr₁.iRs]; hr₂.ndegen[hr₂.iRs]]
    return unique(HRFull(is, js, Rs, vals, numorb(hr₁) + numorb(hr₂); ndegen, inbound=true, dofilter=false))
end
function Base.merge(hr₁::HRFull{T₁}, hr₂::HRFull{T₂}) where {T₁,T₂}
    numorb(hr₁) == numorb(hr₂) || error("Mismatched hr.")
    is = [hr₁.is; hr₂.is]
    js = [hr₁.js; hr₂.js]
    Rs = [hr₁.R[hr₁.iRs]; hr₂.R[hr₂.iRs]]
    vals = [hr₁.vals; hr₂.vals]
    ndegen = [hr₁.ndegen[hr₁.iRs]; hr₂.ndegen[hr₂.iRs]]
    return unique(HRFull(is, js, Rs, vals, numorb(hr₁); ndegen, inbound=true, dofilter=false))
end
function HRFull(path::AbstractMatrix{<:Integer}, vals::AbstractVector{T};
    n::Union{Nothing,Integer}=nothing,
    ndegen::Union{Nothing,AbstractVector{<:Integer}}=nothing,
    inbound::Bool=false, dofilter::Bool=true,
)::HRFull{T} where {T<:Number}

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
        isempty(vals) && return empty(HRFull{T}, 0)
        n = maximum(view(path, 4:5, :))
    else
        isempty(vals) && return empty(HRFull{T}, n)
        inbound || @assert all(i->1 ≤ i ≤ n, view(path, 4:5, :)) "Invalid is detected: i ∈ [1, $n]."
    end

    is = path[4, :]
    js = path[5, :]
    Rs = map(p->ReducedCoordinates{Int}(p[1], p[2], p[3]), eachcol(path))

    return HRFull(is, js, Rs, vals, n; ndegen, inbound=true, dofilter)
end
function HRFull(
    is::AbstractVector{<:Integer},
    js::AbstractVector{<:Integer},
    Rs::AbstractVector{<:AbstractVector{<:Integer}},
    vals::AbstractVector{T},
    n::Integer;
    ndegen::Union{Nothing,AbstractVector{<:Integer}}=nothing,
    inbound::Bool=false, dofilter::Bool=true,
)::HRFull{T} where {T<:Number}

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
    isempty(vals) && return empty(HRFull{T}, n)

    Rs = collect(ReducedCoordinates{Int}, Rs)
    if dofilter
    end
    R = unique(Rs)
    sort!(R)

    R2idx = Dict{ReducedCoordinates{Int},Int}(RR=>i for (i, RR) in enumerate(R))
    iRs = Vector{Int}(undef, length(Rs))
    @inbounds for (i, RR) in enumerate(Rs)
        iRs[i] = R2idx[RR]
    end

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
    end
    return HRFull{T}(n, R, ndegen_uniqueR, is, js, iRs, vals, R2idx)
end
function HRFull(
    is::AbstractVector{<:Integer},
    js::AbstractVector{<:Integer},
    iRs::AbstractVector{<:Integer},
    vals::AbstractVector{T},
    n::Integer,
    R::AbstractVector{<:AbstractVector{<:Integer}};
    ndegen::AbstractVector{<:Integer}=ones(Int, length(R)),
    inbound::Bool=false, dofilter::Bool=true,
)::HRFull{T} where {T<:Number}

    if !inbound
        @assert all(i->1 ≤ i ≤ n, is) "Invalid is detected: i ∈ [1, $n]."
        @assert all(j->1 ≤ j ≤ n, js) "Invalid js detected: j ∈ [1, $n]."
        nR = length(R)
        @assert all(iR->1 ≤ iR ≤ nR, iRs) "Invalid ips detected: iR ∈ [1, $nR]."
        @assert length(ndegen) == nR "Mismatched length of ndegen."
        @assert all(>(0), ndegen) "ndegen should be positive."
        @assert length(is) == length(js) == length(iRs) == length(vals) "Mismatched length of data."
    end
    isempty(vals) && return empty(HRFull{T}, n)

    R = collect(ReducedCoordinates{Int}, R)
    ndegen = collect(Int, ndegen)

    if dofilter
        # filter R depending on if it exists in iRs, then update iRs and ndegen.
        iRs, R, ndegen = _HRFull_filter_R!(iRs, R, ndegen)
    end

    R2idx = Dict{ReducedCoordinates{Int},Int}(RR=>i for (i, RR) in enumerate(R))
    return HRFull{T}(n, R, ndegen, is, js, iRs, vals, R2idx)
end
function _HRFull_filter_R!(iRs, R, ndegen)
    R_exist = unique(iRs)
    sort!(R_exist)
    new_R = unique(view(R, R_exist))
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
"""
	spinorize(hr::HRFull{T}; mode = conj) -> HRFull{T}

Generate a spinful hr from a spinless hr.
Acctually, this function just makes a copy by `mode`, it won't introduce new hopping.
"""
function spinorize(hr::HRFull{T}; mode=conj) where {T}
    is = [hr.is; hr.is .+ numorb(hr)]
    js = [hr.js; hr.js .+ numorb(hr)]
    iRs = [hr.iRs; hr.iRs]
    vals = [hr.vals; mode.(hr.vals)]
    return HRFull{T}(numorb(hr) * 2, hr.R, hr.ndegen, is, js, iRs, vals, hr.R2idx)
end
"""
	translate(hr::HRFull{T}, centres::Pair{<:Integer, <:AbstractVector{<:Integer}}...) -> HRFull{T}

You only need to input the wannier centres R which isn't in the [0, 0, 0] cell.

See also [`translate!`](@ref).
"""
function translate(hr::HRFull{T}, aim_centres::Pair{<:Integer,<:AbstractVector{<:Integer}}...) where {T}
    # work for infinite
    centres = [ReducedCoordinates{Int}(0, 0, 0) for _ in 1:numorb(hr)]
    for (iw, R) in aim_centres
        centres[iw] = ReducedCoordinates{Int}(R)
    end
    Rs = Vector{ReducedCoordinates{Int}}(undef, length(hr))
    vals = Vector{T}(undef, length(hr))
    for ih in eachindex(Rs)
        i, j, iR = hr.is[ih], hr.js[ih], hr.iRs[ih]
        Rs[ih] = hr.R[iR] + centres[j] - centres[i]
        # vals[ih] = hr.vals[ih] / hr.ndegen[iR]
    end
    return HRFull(copy(hr.is), copy(hr.js), Rs, vals, numorb(hr); inbound=true, dofilter=false)
end
function translate(hr::HRFull{T}, lattice::Union{AbstractMatrix{<:Real},Lattice{<:Real}}, grid,
    aim_centres::Pair{<:Integer,<:AbstractVector{<:Integer}}...;
    atol::Real=1e-6, rtol::Real=atol>0 ? 0 : √eps()) where {T}
    # work for finite
    # First, obtain independent values!
    grid_range = gridindex(Val(:range), grid)
    Rref = [ReducedCoordinates(i, j, k) for i ∈ grid_range[1], j ∈ grid_range[2], k ∈ grid_range[3]]
    Rref = reshape(Rref, :)
    Rrefs2idx = Dict{ReducedCoordinates{Int},Int}(RR=>iR for (iR, RR) in enumerate(Rref))
    # orbital centres
    norb = numorb(hr)
    centres = [ReducedCoordinates{Int}(0, 0, 0) for _ in 1:norb]
    for (iw, R) in aim_centres
        centres[iw] = ReducedCoordinates{Int}(R)
    end

    # if linear storage, then need dict
    # Repeated elements will be overwritten.
    vals_array = zeros(T, norb, norb, length(Rref))
    for (i, j, iR, v) in enumerate(hr)
        newR = hr.R[iR] + centres[j] - centres[i]
        newRref = ReducedCoordinates{Int}(map(mod, newR, grid_range))
        iRref = Rrefs2idx[newRref]
        vals_array[i, j, iRref] = v
    end

    # obtain ndegen
    R_kinds = map(RR -> _grid_degeneracy(RR, grid, lattice; atol, rtol), Rref)
    ndegen_kinds = length.(R_kinds)
    nR = sum(ndegen_kinds)
    allR = Vector{ReducedCoordinates{Int}}(undef, 0)
    sizehint!(allR, nR)
    ndegen = Vector{Int}(undef, 0)
    sizehint!(ndegen, nR)
    for R_kind in R_kinds
        ndegen_kind = length(R_kind)
        append!(allR, R_kind)
        append!(ndegen, fill(ndegen_kind, ndegen_kind))
    end
    perm = sortperm(allR)
    permute!(allR, perm)
    permute!(ndegen, perm)
    allR2idx = Dict{ReducedCoordinates{Int},Int}(RR=>iR for (iR, RR) in enumerate(allR))

    # construct data.
    nR = length(allR)
    N = norb * norb * nR
    is = Vector{Int}(undef, 0)
    js = Vector{Int}(undef, 0)
    iRs = Vector{Int}(undef, 0)
    vals = Vector{T}(undef, 0)
    sizehint!(is, N)
    sizehint!(js, N)
    sizehint!(iRs, N)
    sizehint!(vals, N)
    for iR in eachindex(Rref)
        iallRs = map(RR->allR2idx[RR], R_kinds[iR])
        for j in Base.OneTo(norb), i in Base.OneTo(norb)
            v = vals_array[i, j, iR]
            iszero(v) && continue
            for iallR in iallRs
                push!(is, i)
                push!(js, j)
                push!(iRs, iallR)
                push!(vals, v)
            end
        end
    end

    hr = HRFull(is, js, iRs, vals, norb, allR; ndegen, inbound=true, dofilter=true)
    return sort!(hr)
end
"""
	translate!(hr::HRFull{T}, centres::Pair{<:Integer, <:AbstractVector{<:Integer}}...) -> HRFull{T}

You only need to input the wannier centres R which isn't in the [0, 0, 0] cell.

See also [`translate`](@ref).
"""
function translate!(hr::HRFull{T}, aim_centres::Pair{<:Integer,<:AbstractVector{<:Integer}}...) where {T}
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
    empty!(hr.ndegen)
    append!(hr.R, unique(Rs))
    sort!(hr.R)
    for (iR, RR) in enumerate(hr.R)
        hr.R2idx[RR] = iR
    end
    append!(hr.ndegen, fill(1, hr.nR))
    for (ih, RR) in enumerate(Rs)
        hr.iRs[ih] = hr.R2idx[RR]
    end
    return hr
end
