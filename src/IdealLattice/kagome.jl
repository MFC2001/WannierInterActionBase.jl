function IdealCell(::Val{:kagome}; a=1, b=a, c=20)

    l = a + b
    𝐚 = Vec3(l, 0, 0)
    𝐛 = l * Vec3(-1 / 2, √3 / 2, 0)
    𝐜 = Vec3(0, 0, c)

    lattice = Lattice([𝐚 𝐛 𝐜])

    location2 = [b, 0, c / 2]
    location1 = [b - a / 2, a * √3 / 2, c / 2]
    location3 = [b + a / 2, a * √3 / 2, c / 2]

    location = CartesianCoordinates.([location1, location2, location3])
    location = map(x -> lattice \ x, location)

    return Cell(lattice, location, "ReducedCoordinates"; name=["A", "B", "C"], period=Bool[1, 1, 0])
end
function IdealHR(::Val{:kagome}; t₀=0,
    t₁=1, t₁_intra=t₁, t₁_inter=t₁,
    θ₁ₕ=0, θ₁ₕ_intra=θ₁ₕ, θ₁ₕ_inter=θ₁ₕ)
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
    return HR(allpath, allvalue; n=3)
end
