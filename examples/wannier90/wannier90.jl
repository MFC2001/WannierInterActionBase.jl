
using LinearAlgebra
using WannierInterActionBase

a₁ = [1, 0, 0]
a₂ = [0, 1, 0]
a₃ = [0, 0, 1]
lattice = Lattice(a₁, a₂, a₃)
cell = Cell(lattice, [[0, 0, 0]]; name=["X"])

mp_grid = [7, 7, 1]
kgrid = RedKgrid(MonkhorstPack(mp_grid))
write("./wannier90.win", cell, wannier90_win;
    num_iter=2000,
    num_bands=3,
    dis_win_min=0,
    dis_win_max=20,
    dis_froz_min=0,
    dis_froz_max=20,
    num_wann=3,
    mp_grid,
    kpoints=kgrid,
)

(nnkpts, weights) = kgrid_shell_automatic(kgrid, reciprocal(lattice))
nkb = size(nnkpts, 2)

band = map(k->eigen(H(k)), kgrid)

# band can be a Vector{Eigen} or a Matrix{Real}(nb, nk)
write("./wannier90.eig", band, wannier90_eig)
# M isa Array{Number}(nb, nb, nkb)
write("./wannier90.mmn", M, wannier90_mmn; nnkpts)
# A isa Array{Number}(nb, nw, nk)
write("./wannier90.amn", A, wannier90_amn)

