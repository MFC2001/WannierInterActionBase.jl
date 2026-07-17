
using WannierInterActionBase

function Gauss_αrcut(x, α; ϵ = 1, atol = 1e-8)
	V = RealCoulomb(:point; ϵ)
	VG = RealCoulomb(:gauss; ϵ, α)
	i = findfirst(x) do xx
		isapprox(VG(xx), V(xx); atol)
	end
	rcut = x[i]
	αrcut = rcut * α
	return V(rcut), rcut, αrcut
end

Gauss_αrcut(0.001:0.001:100, 3)
Gauss_αrcut(0.0001:0.0001:100, 15) # rcut = 0.3, αrcut = 4.5
