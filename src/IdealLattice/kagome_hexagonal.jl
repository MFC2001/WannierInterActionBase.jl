function IdealCell(::Val{:kagome_hexagonal}; a=1, b=a, c=20)

    l = a + b
    𝐚 = Vec3(l, 0, 0)
    𝐛 = l * Vec3(-1 / 2, √3 / 2, 0)
    𝐜 = Vec3(0, 0, c)

    lattice = Lattice([𝐚 𝐛 𝐜])

    location2 = [b, 0, c / 2]
    location1 = [b - a / 2, a * √3 / 2, c / 2]
    location3 = [b + a / 2, a * √3 / 2, c / 2]

    location4 = [b, a * √3 / 3, c / 2]
    location5 = [b / 2 - a / 2, a * √3 / 2 + b * √3 / 6, c / 2]

    location = CartesianCoordinates.([location1, location2, location3, location4, location5])
    location = map(x -> lattice \ x, location)

    return Cell(lattice, location, "ReducedCoordinates"; name=["A", "B", "C", "a", "b"], period=Bool[1, 1, 0])
end
function IdealHR(::Val{:kagome_hexagonal};
    t₀=0, t₀_a=t₀, t₀_b=t₀,
    t₁=1, t₁_intra=t₁, t₁_inter=t₁,
    t₁_a=0, t₁_b=0,
    t₁_ab=0, t₂_ab=0,
    θ₁ₕ=0, θ₁ₕ_intra=θ₁ₕ, θ₁ₕ_inter=θ₁ₕ,
    θ₂ₕ=0)
    allpath = Matrix{Int}(undef, 5, 0)
    allvalue = Vector{Float64}(undef, 0)
    if !iszero(t₀)
        path_t₀ = [
            0 0 0 1 1;
            0 0 0 2 2;
            0 0 0 3 3;
        ]
        value_t₀ = fill(t₀, size(path_t₀, 1))
        allpath = [allpath transpose(path_t₀)]
        allvalue = [allvalue; value_t₀]
    end
    if !iszero(t₀_a)
        path_t₀_a = [0 0 0 4 4]
        allpath = [allpath transpose(path_t₀_a)]
        allvalue = [allvalue; t₀_a]
    end
    if !iszero(t₀_b)
        path_t₀_b = [0 0 0 5 5]
        allpath = [allpath transpose(path_t₀_b)]
        allvalue = [allvalue; t₀_b]
    end
    if !iszero(t₁_intra)
        path_t₁_intra = [
            0 0 0 1 2 -1;
            0 0 0 1 3 1;
            0 0 0 2 1 1;
            0 0 0 2 3 -1;
            0 0 0 3 1 -1;
            0 0 0 3 2 1;
        ]
        if iszero(θ₁ₕ_intra)
            phase_t₁_intra = ones(Int, size(path_t₁_intra, 1))
        else
            phase_t₁_intra = cis.(θ₁ₕ_intra .* path_t₁_intra[:, 6])
        end
        path_t₁_intra = path_t₁_intra[:, 1:5]
        value_t₁_intra = t₁_intra .* phase_t₁_intra
        allpath = [allpath transpose(path_t₁_intra)]
        allvalue = [allvalue; value_t₁_intra]
    end
    if !iszero(t₁_inter)
        path_t₁_inter = [
            0 -1 0 1 2 -1;
            -1 -1 0 1 3 1;
            0 1 0 2 1 1;
            -1 0 0 2 3 -1;
            1 1 0 3 1 -1;
            1 0 0 3 2 1;
        ]
        if iszero(θ₁ₕ_inter)
            phase_t₁_inter = ones(Int, size(path_t₁_inter, 1))
        else
            phase_t₁_inter = cis.(θ₁ₕ_inter .* path_t₁_inter[:, 6])
        end
        path_t₁_inter = path_t₁_inter[:, 1:5]
        value_t₁_inter = t₁_inter .* phase_t₁_inter
        allpath = [allpath transpose(path_t₁_inter)]
        allvalue = [allvalue; value_t₁_inter]
    end
    if !iszero(t₁_a)
        path_t₁_a = [
            0 0 0 4 1;
            0 0 0 4 2;
            0 0 0 4 3;
            0 0 0 1 4;
            0 0 0 2 4;
            0 0 0 3 4;
        ]
        value_t₁_a = fill(t₁_a, size(path_t₁_a, 1))
        allpath = [allpath transpose(path_t₁_a)]
        allvalue = [allvalue; value_t₁_a]
    end
    if !iszero(t₁_b)
        path_t₁_b = [
            0 1 0 5 1;
            0 0 0 5 2;
            -1 0 0 5 3;
            0 -1 0 1 5;
            0 0 0 2 5;
            1 0 0 3 5
        ]
        value_t₁_b = fill(t₁_b, size(path_t₁_b, 1))
        allpath = [allpath transpose(path_t₁_b)]
        allvalue = [allvalue; value_t₁_b]
    end
    if !iszero(t₁_ab)
        path_t₁_ab = [
            0 0 0 4 5;
            1 0 0 4 5;
            0 -1 0 4 5;
            0 0 0 5 4;
            -1 0 0 5 4;
            0 1 0 5 4;
        ]
        value_t₁_ab = fill(t₁_ab, size(path_t₁_ab, 1))
        allpath = [allpath transpose(path_t₁_ab)]
        allvalue = [allvalue; value_t₁_ab]
    end
    if !iszero(t₂_ab)
        path_t₂_ab = [
            1 1 0 4 4 1;
            1 0 0 4 4 -1;
            0 1 0 4 4 -1;
            -1 0 0 4 4 1;
            0 -1 0 4 4 1;
            -1 -1 0 4 4 -1;
            1 1 0 5 5 -1;
            1 0 0 5 5 1;
            0 1 0 5 5 1;
            -1 0 0 5 5 -1;
            0 -1 0 5 5 -1;
            -1 -1 0 5 5 1;
        ]
        if iszero(θ₂ₕ)
            phase_t₂_ab = ones(Int, size(path_t₂_ab, 1))
        else
            phase_t₂_ab = cis.(θ₂ₕ .* path_t₂_ab[:, 6])
        end
        path_t₂_ab = path_t₂_ab[:, 1:5]
        value_t₂_ab = t₂_ab .* phase_t₂_ab
        allpath = [allpath transpose(path_t₂_ab)]
        allvalue = [allvalue; value_t₂_ab]
    end
    return HR(allpath, allvalue; indices=collect(1:5), dosort=true)
end
