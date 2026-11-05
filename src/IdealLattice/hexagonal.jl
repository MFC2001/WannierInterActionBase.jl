function IdealCell(::Val{:hexagonal}; a=1, c=20)

    l = √3 * a
    𝐚 = [l, 0, 0]
    𝐛 = l * [-1 / 2, √3 / 2, 0]
    𝐜 = [0, 0, c]

    lattice = Lattice([𝐚 𝐛 𝐜])

    location1 = [0, a, c / 2]
    location2 = [l / 2, a / 2, c / 2]

    location = CartesianCoordinates.([location1, location2])
    location = map(x -> lattice \ x, location)

    return Cell(lattice, location, "Reduced"; name=["A", "B"], period=Bool[1, 1, 0])
end
function IdealHR(::Val{:hexagonal}; t₀=0, tAB=0, t₀_A=t₀+tAB, t₀_B=t₀-tAB,
    t₁=1, t₂=0, θ₂ₕ=0)
    allpath = Matrix{Int}(undef, 5, 0)
    allvalue = Vector{Float64}(undef, 0)
    if !iszero(t₀_A)
        path_t₀_A = [0 0 0 1 1]
        allpath = [allpath transpose(path_t₀_A)]
        allvalue = [allvalue; t₀_A]
    end
    if !iszero(t₀_B)
        path_t₀_B = [0 0 0 2 2]
        allpath = [allpath transpose(path_t₀_B)]
        allvalue = [allvalue; t₀_B]
    end
    if !iszero(t₁)
        path_t₁ = [
            0 0 0 1 2;
            -1 0 0 1 2;
            0 1 0 1 2;
            0 0 0 2 1;
            1 0 0 2 1;
            0 -1 0 2 1
        ]
        value_t₁ = fill(t₁, size(path_t₁, 1))
        allpath = [allpath transpose(path_t₁)]
        allvalue = [allvalue; value_t₁]
    end
    if !iszero(t₂)
        path_t₂ = [
            1 1 0 1 1 -1;
            1 0 0 1 1 1;
            0 1 0 1 1 1;
            -1 0 0 1 1 -1;
            0 -1 0 1 1 -1;
            -1 -1 0 1 1 1;
            1 1 0 2 2 1;
            1 0 0 2 2 -1;
            0 1 0 2 2 -1;
            -1 0 0 2 2 1;
            0 -1 0 2 2 1;
            -1 -1 0 2 2 -1
        ]
        if iszero(θ₂ₕ)
            phase_t₂ = ones(Int, size(path_t₂, 1))
        else
            phase_t₂ = cis.(θ₂ₕ .* path_t₂[:, 6])
        end
        path_t₂ = path_t₂[:, 1:5]
        value_t₂ = t₂ .* phase_t₂
        allpath = [allpath transpose(path_t₂)]
        allvalue = [allvalue; value_t₂]
    end
    return HR(allpath, allvalue; indices=collect(1:2), dosort=true)
end
