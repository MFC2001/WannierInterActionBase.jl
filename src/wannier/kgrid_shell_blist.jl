function kgrid_shell_blist(kgrid, lattice; kwargs...)
    lattice = Lattice(parent(lattice))
    rlattice = reciprocal(lattice)
    return kgrid_shell_blist(kgrid, rlattice; kwargs...)
end
function kgrid_shell_blist(kgrid::MonkhorstPack, rlattice::ReciprocalLattice; kwargs...)
    kgrid = RedKgrid(kgrid)
    return kgrid_shell_blist(kgrid, rlattice; kwargs...)
end
function kgrid_shell_blist(kgrid::RedKgrid, rlattice::ReciprocalLattice;
    nsupercell::Integer=2, shell_tol::Real=1e-6, orth_tol::Real=1e-6)

    @assert shell_tol > 0 "`shell_tol` should be positive!"
    @assert orth_tol > 0 "`orth_tol` should be positive!"

    I = _wannier_get_nsupercell(rlattice, nsupercell)

    Nk = length(kgrid)
    kgrid_car = map(k -> rlattice * k, kgrid)

    #计算并保存所有b
    Nb = Nk * (2 * I[1] + 1) * (2 * I[2] + 1) * (2 * I[3] + 1)
    b_all = Vector{Vec3{Float64}}(undef, Nb)
    bindex = Matrix{Int}(undef, 4, Nb)
    nb = 0
    for n in (-I[3]):I[3], m in (-I[2]):I[2], l in (-I[1]):I[1]
        lmn = rlattice * [l, m, n]
        for (ki, kcar) in enumerate(kgrid_car)
            nb += 1
            b_all[nb] = Vec3(kcar + lmn)
            bindex[:, nb] = [ki, l, m, n]
        end
    end

    #search shell.
    dnn_all = map(b -> norm(b), b_all)
    dnn_index = sortperm(dnn_all)
    dnn_index = dnn_index[2:end]
    dnn_all = dnn_all[dnn_index]

    dnn = Vector{Float64}(undef, 1)
    shell_bindex = Vector{Vector{Int}}(undef, 1)
    dnn[1] = dnn_all[1]
    shell_bindex[1] = Int[]
    for (dnn_i, dnn_v) in zip(dnn_index, dnn_all)
        if isapprox(dnn_v, dnn[end]; atol=shell_tol)
            push!(shell_bindex[end], dnn_i)
        elseif dnn_v > dnn[end] + shell_tol
            push!(dnn, dnn_v)
            push!(shell_bindex, [dnn_i])
        else
            error("Shell search failed: try adjusting `shell_tol` or fixing basis symmetry.")
        end
    end
    multi = length.(shell_bindex)

    #calculate weights of shells
    aimshell_index = Vector{Int}(undef, 0)
    aimb_index = Vector{Int}(undef, 0)
    Amat = Matrix{Float64}(undef, 6, 0)
    issatisfy = false
    for ishell in 1:length(dnn)
        _wannier_check_shell_parallel(b_all, aimb_index, shell_bindex[ishell]) && continue
        T = sum(b_all[shell_bindex[ishell]]) do b
            [b[1]^2, b[2]^2, b[3]^2, b[1] * b[2], b[2] * b[3], b[3] * b[1]]
        end
        F = svd([Amat T])
        any(abs.(F.S) .< 1e-5) && continue
        Amat = [Amat T]
        push!(aimshell_index, ishell)
        append!(aimb_index, shell_bindex[ishell])
        weights = F.V * Diagonal(1 ./ F.S) * F.U' * [1, 1, 1, 0, 0, 0]
        if isapprox(Amat * weights, [1, 1, 1, 0, 0, 0]; atol=orth_tol)
            issatisfy = true
            break
        end
    end

    if !issatisfy
        error("Can't satisfy weighted orthogonality, try to increase `nsupercell` or decrease `shell_tol`.")
    end

    blist = map(aimb_index) do ib
        ReducedCoordinates(kgrid[bindex[1, ib]] + bindex[2:4, ib])
    end

    F = svd(Amat)
    weights_shell = F.V * Diagonal(1 ./ F.S) * F.U' * [1, 1, 1, 0, 0, 0]
    aim_multi = multi[aimshell_index]
    weights = Vector{Float64}(undef, 0)
    for (i, w) in enumerate(weights_shell)
        append!(weights, fill(w, aim_multi[i]))
    end

    return blist, weights
end
function _wannier_check_shell_parallel(b_all, aimb_index, newb_index)
    for newb_i in newb_index
        newb = b_all[newb_i]
        norm_newb = norm(newb)
        for aimb_i in aimb_index
            cosθ = (newb ⋅ b_all[aimb_i]) / (norm_newb * norm(b_all[aimb_i]))
            if cosθ ≈ 1
                return true
            end
        end
    end
    return false
end
function _wannier_get_nsupercell(rlattice, nsupercell)

    a₁ = rlattice[:, 1]
    a₂ = rlattice[:, 2]
    a₃ = rlattice[:, 3]

    V = ((a₁ × a₂) ⋅ a₃)

    h₁ = V / norm(a₂ × a₃)
    h₂ = V / norm(a₃ × a₁)
    h₃ = V / norm(a₁ × a₂)

    hmax = max(h₁, h₂, h₃)
    I = ceil.((nsupercell * hmax) ./ [h₁, h₂, h₃])

    return Int.(I)
end
