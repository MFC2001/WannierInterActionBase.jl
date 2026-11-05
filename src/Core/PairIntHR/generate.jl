export generate_PairIntHR
function generate_PairIntHR(pair::AbstractVector{<:WanIntPairVecType},
    lattice, orblocat, grid; pv, V, atol::Real=1e-6)

    lattice = Lattice(parent(lattice))
    orblocat = ReducedCoordinates.(orblocat)


    pair_r₀ = map(pair) do p
        lattice * (orblocat[p[2]] + p[3] + orblocat[p[1]]) / 2
    end
    pair_Δr = map(pair) do p
        Δr = lattice * (orblocat[p[2]] + p[3] - orblocat[p[1]])
        norm(Δr)
    end

    pair_ρ = pv.(pair_Δr)

    grid = gridindex(grid)
    grid_car = map(R->lattice * R, grid)

    np = length(pair)
    value = [pair_ρ[p1] * pair_ρ[p2] * V(pair_r₀[p2] + R - pair_r₀[p1]) for p1 in 1:np, p2 in 1:np, R in grid_car]

    return PairIntHR(value, pair, grid; nw=length(orblocat), atol, dosort=true, dofilter=true)
end