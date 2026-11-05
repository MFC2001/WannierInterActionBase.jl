gridindex(sym::Symbol, args...; kwargs...) = gridindex(Val(sym), args...; kwargs...)
function gridindex(::Val{:range}, grid)
    start = map(x->-x, floor.(Int, (grid .- 1) .// 2))
    stop = ceil.(Int, (grid .- 1) .// 2)
    return ntuple(i->start[i]:stop[i], Val(length(grid)))
end
"""
	gridindex(grid) -> Vector{ReducedCoordinates{Int}}

The examples of `grid` are [3,3,1], [-1:1, 0:0, -1:1], (3,3,3) or (-1:1,-1:1,-1:1).
Return a vector like [R[0,0,0], R[1,2,3],...].
The interval is open on the left and closed on the right, e.g.(-grid/2, grid/2].
"""
function gridindex(grid::AbstractVector{<:Integer})::Vector{ReducedCoordinates{Int}}
    length(grid) == 3 || throw(ArgumentError("should input 3 integers vector for grid"))
    grid_range = gridindex(Val(:range), grid)
    return vec([ReducedCoordinates(i, j, k) for i ∈ grid_range[1], j ∈ grid_range[2], k ∈ grid_range[3]])
end
function gridindex(grid::NTuple{3,<:Integer})::Vector{ReducedCoordinates{Int}}
    grid_range = gridindex(Val(:range), grid)
    return vec([ReducedCoordinates(i, j, k) for i ∈ grid_range[1], j ∈ grid_range[2], k ∈ grid_range[3]])
end
function gridindex(grid::AbstractVector{<:UnitRange{<:Integer}})::Vector{ReducedCoordinates{Int}}
    length(grid) == 3 || throw(ArgumentError("should input 3 integers vector for grid"))
    all(x -> step(x) == 1, grid) || throw(ArgumentError("step should be 1 for each dimension in grid"))
    grid = collect(UnitRange{Int}, grid)
    return vec([ReducedCoordinates(i, j, k) for i ∈ grid[1], j ∈ grid[2], k ∈ grid[3]])
end
function gridindex(grid::NTuple{3,<:UnitRange{<:Integer}})::Vector{ReducedCoordinates{Int}}
    all(x -> step(x) == 1, grid) || throw(ArgumentError("step should be 1 for each dimension in grid"))
    grid = UnitRange{Int}.(grid)
    return vec([ReducedCoordinates(i, j, k) for i ∈ grid[1], j ∈ grid[2], k ∈ grid[3]])
end
"""
	gridindex(lattice, maxlattdist::Real, period::AbstractVector{Bool}) -> Vector{ReducedCoordinates{Int}}

Get all the lattice vectors whose norms are smaller than maxlattdist.
"""
function gridindex(lattice::Union{AbstractMatrix{<:Real},AbstractLattice{<:Real}}, maxdist::Real, period::AbstractVector{Bool})::Vector{ReducedCoordinates{Int}}

    (N1, N2, N3) = fractional_bounds(lattice, maxdist; atol=1e-4)

    period = Vec3{Bool}(period)
    p = count(period)
    if p == 3
        R = gridindex(((-N1):N1, (-N2):N2, (-N3):N3))
    elseif p == 2
        if !period[1]
            R = gridindex((0:0, (-N2):N2, (-N3):N3))
        elseif !period[2]
            R = gridindex(((-N1):N1, 0:0, (-N3):N3))
        else
            R = gridindex(((-N1):N1, (-N2):N2, 0:0))
        end
    elseif p == 1
        if period[1]
            R = gridindex(((-N1):N1, 0:0, 0:0))
        elseif period[2]
            R = gridindex((0:0, (-N2):N2, 0:0))
        else
            R = gridindex((0:0, 0:0, (-N3):N3))
        end
    else
        R = [ReducedCoordinates(0, 0, 0)]
    end

    ldot = parent(lattice)
    ldot = Mat3(transpose(ldot) * ldot)
    maxdist2 = maxdist^2
    R = filter(RR -> RR ⋅ (ldot * RR) < maxdist2, R)

    return R
end
function gridindex(::Val{:WignerSeitz}, lattice::Union{AbstractMatrix{<:Real},AbstractLattice{<:Real}}, grid;
    atol::Real=1e-6, rtol::Real=atol>0 ? 0 : √eps())
    Rrefs = gridindex(grid)
    grid = Vec3{Int}(grid)
    R_kinds = map(RR -> _grid_degeneracy(RR, grid, lattice; atol, rtol), Rrefs)
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
    return allR, ndegen
end
function _grid_degeneracy(Rref, grid, lattice; atol::Real=1e-6, rtol::Real=atol>0 ? 0 : √eps())
    Rref_new = Vector{Int}(undef, 3)
    for i in 1:3
        x = mod(Rref[i], grid[i])
        if 2 * x > grid[i]
            Rref_new[i] = x - grid[i]
        else
            Rref_new[i] = x
        end
    end
    Rref = ReducedCoordinates{Int}(Rref_new)

    if iszero(Rref)
        return [Rref]
    end

    lattice = Lattice(parent(lattice))
    Rref_car = lattice * Rref
    Rref_norm = norm(Rref_car)

    lattice_grid = Matrix{Float64}(undef, 3, 3)
    lattice_grid[:, 1] = lattice[:, 1] .* grid[1]
    lattice_grid[:, 2] = lattice[:, 2] .* grid[2]
    lattice_grid[:, 3] = lattice[:, 3] .* grid[3]
    lattice_grid = Lattice(lattice_grid)

    # 这里的atol只是避免结果为整数的情况，宽松一点只会导致增加一点点计算量而已。
    (N1, N2, N3) = fractional_bounds(lattice_grid, 2.0 * Rref_norm; atol=1e-4)
    N = (2*N1+1)*(2*N2+1)*(2*N3+1)
    if N == 1
        return [Rref]
    end

    ldot = parent(lattice)
    ldot = Mat3(transpose(ldot) * ldot)
    min_d2 = Inf
    allR = Vector{ReducedCoordinates{Int}}(undef, N)
    allR_norm2 = Vector{Float64}(undef, N)
    nn = 0
    for n3 in (-N3):N3, n2 in (-N2):N2, n1 in (-N1):N1
        Rnew = Rref + ReducedCoordinates{Int}(n1*grid[1], n2*grid[2], n3*grid[3])
        Rnew_norm2 = Rnew ⋅ (ldot * Rnew)
        if Rnew_norm2 < min_d2
            min_d2 = Rnew_norm2
        end
        nn+=1
        allR[nn] = Rnew
        allR_norm2[nn] = Rnew_norm2
    end

    tol = max(atol, rtol * min_d2)
    shortest = allR[allR_norm2 .≤ (min_d2+tol)]

    sort!(shortest)
    return shortest
end
# function gridindex(lattice, maxlattdist, ::Val{3})

#     𝐚 = lattice[begin:(begin+2)]
#     𝐛 = lattice[(begin+3):(begin+5)]
#     𝐜 = lattice[(begin+6):end]

#     Ω = abs((𝐚 × 𝐛) ⋅ 𝐜)

#     h₁ = Ω / norm(𝐛 × 𝐜)
#     h₂ = Ω / norm(𝐜 × 𝐚)
#     h₃ = Ω / norm(𝐚 × 𝐛)

#     grid = gridindex(Int.(cld.(maxlattdist * 1.2, [h₁, h₂, h₃])) * 2 .+ 1)
#     maxlattdist2 = maxlattdist^2
#     ldot = parent(lattice)
#     ldot = Mat3(transpose(ldot) * ldot)
#     grid = filter(R -> R ⋅ (ldot * R) < maxlattdist2, grid)

#     return grid
# end
# function gridindex(lattice, maxlattdist, xyindex, ::Val{2})

#     𝐚 = lattice[begin:(begin+2)]
#     𝐛 = lattice[(begin+3):(begin+5)]
#     𝐜 = lattice[(begin+6):end]

#     Ω = abs((𝐚 × 𝐛) ⋅ 𝐜)

#     h₁ = Ω / norm(𝐛 × 𝐜)
#     h₂ = Ω / norm(𝐜 × 𝐚)
#     h₃ = Ω / norm(𝐚 × 𝐛)


#     grid_3D = Int.(cld.(maxlattdist * 1.2, [h₁, h₂, h₃])) * 2 .+ 1
#     grid = [1, 1, 1]
#     grid[xyindex] .= grid_3D[xyindex]

#     grid = gridindex(grid)
#     maxlattdist2 = maxlattdist^2
#     ldot = parent(lattice)
#     ldot = Mat3(transpose(ldot) * ldot)
#     grid = filter(R -> R ⋅ (ldot * R) < maxlattdist2, grid)

#     return grid
# end
# function gridindex(lattice, maxlattdist, zindex, ::Val{1})

#     𝐚 = lattice[begin:(begin+2)]
#     𝐛 = lattice[(begin+3):(begin+5)]
#     𝐜 = lattice[(begin+6):end]

#     Ω = abs((𝐚 × 𝐛) ⋅ 𝐜)

#     h₁ = Ω / norm(𝐛 × 𝐜)
#     h₂ = Ω / norm(𝐜 × 𝐚)
#     h₃ = Ω / norm(𝐚 × 𝐛)

#     grid_3D = Int.(cld.(maxlattdist * 1.2, [h₁, h₂, h₃])) * 2 .+ 1
#     grid = [1, 1, 1]
#     grid[zindex] = grid_3D[zindex]

#     grid = gridindex(grid)
#     maxlattdist2 = maxlattdist^2
#     ldot = parent(lattice)
#     ldot = Mat3(transpose(ldot) * ldot)
#     grid = filter(R -> R ⋅ (ldot * R) < maxlattdist2, grid)

#     return grid
# end

