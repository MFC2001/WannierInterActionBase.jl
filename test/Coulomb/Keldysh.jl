
using WannierInterActionBase

function Keldysh_r₀(x, r₀; atol = 1e-8)
	V = RealCoulomb(:point)
	VK = RealCoulomb(:keldysh; r₀)
	i = findfirst(x) do xx
		isapprox(VK(xx), V(xx); atol)
	end
	rcut = x[i]
	return V(rcut), rcut, rcut / r₀
end

Keldysh_r₀(1:1:10000, 1; atol = 1e-6)
Keldysh_r₀(1:1:10000, 10)
