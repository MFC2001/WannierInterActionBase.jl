"""
	kgrid_shell_automatic(kgrid, lattice; kwargs...) -> (nnkpts, weights)

Select the shells of neighbouring k-points needed to build a real-space Wannier
Hamiltonian by an exact inverse Fourier transform, and return the neighbour table
`nnkpts` (5 rows holding `ik`, `ikb`, `l`, `m`, `n`) together with the shell weights.

`nsupercell` sets the size of the supercell search, `shell_tol` groups vectors of equal
length and `orth_tol` controls the weighted-orthogonality check; an error is thrown when
no shell combination satisfies that constraint.
"""
function kgrid_shell_automatic(kgrid, lattice; kwargs...)
    lattice = Lattice(parent(lattice))
    rlattice = reciprocal(lattice)
    return kgrid_shell_automatic(kgrid, rlattice; kwargs...)
end
function kgrid_shell_automatic(kgrid::MonkhorstPack, rlattice::ReciprocalLattice; kwargs...)
    kgrid = RedKgrid(kgrid)
    return kgrid_shell_automatic(kgrid, rlattice; kwargs...)
end
function kgrid_shell_automatic(kgrid::RedKgrid, rlattice::ReciprocalLattice;
    nsupercell::Integer=3, shell_tol::Real=1e-6, orth_tol::Real=1e-6)

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

    aim_dnn = dnn[aimshell_index]
    aim_multi = multi[aimshell_index]

    aim_nb = sum(aim_multi)
    nnkpts = Matrix{Int}(undef, 5, Nk * aim_nb)
    weights_all = Vector{Float64}(undef, 0)
    sizehint!(weights_all, Nk * aim_nb)
    Amat = Matrix{Float64}(undef, 6, length(aim_dnn))
    for (ki, kcar) in enumerate(kgrid_car)
        b_all_k = map(b -> b - kcar, b_all)
        dnn_all = map(b -> norm(b), b_all_k)

        b_neighbour = Vector{Int}(undef, 0)
        for (i, d) in enumerate(aim_dnn)
            I = findall(dist -> isapprox(dist, d; atol=shell_tol), dnn_all)
            length(I) == aim_multi[i] || error("Searched shells don't have generalization.")
            append!(b_neighbour, I)
            Amat[:, i] = sum(b_all_k[I]) do b
                [b[1]^2, b[2]^2, b[3]^2, b[1] * b[2], b[2] * b[3], b[3] * b[1]]
            end
        end
        F = svd(Amat)
        weights = F.V * Diagonal(1 ./ F.S) * F.U' * [1, 1, 1, 0, 0, 0]
        isapprox(Amat * weights, [1, 1, 1, 0, 0, 0]; atol=orth_tol) || error("Searched shells don't have generalization.")

        for (i, w) in enumerate(weights)
            append!(weights_all, fill(w, aim_multi[i]))
        end
        nnkpts[:, ((ki-1)*aim_nb+1):(ki*aim_nb)] = [fill(ki, 1, aim_nb); bindex[:, b_neighbour]]
    end

    # 每个ik都有相同个数的b，且连续存储。
    return nnkpts, weights_all
end
