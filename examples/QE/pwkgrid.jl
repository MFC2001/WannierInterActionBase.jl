
using WannierInterActionBase

kpoints = RedKgrid(MonkhorstPack([23, 23, 1]))
#sig2wan.x of BerkeleyGW need k in [0, 1]
# kpoints = map(kpoints.kdirect) do k
# 	mod.(k, 1)
# end
#The txt of kpoints can be copy to .win file.
write("./pw.in", zero(Cell), QEpw; kpoints)
