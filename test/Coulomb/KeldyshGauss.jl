
using WannierInterActionBase

function KeldyshGauss_αrcut(x, α; r₀ = 1, atol = 1e-8)
	VK = RealCoulomb(:keldysh; r₀)
	VKG = RealCoulomb(:keldysh_gauss; r₀, α)
	i = findfirst(x) do xx
		isapprox(VKG(xx), VK(xx); atol)
	end
	rcut = x[i]
	αrcut = rcut * α
	return VK(rcut), rcut, αrcut
end

KeldyshGauss_αrcut(1:1:1000, 1; r₀ = 1, atol = 1e-6)



# using SpecialFunctions
# using QuadGK
# function fint(q, α)
# 	integrand(kz) = exp(-kz^2/(4 * α^2))/(kz^2 + q^2)
# 	result, err = quadgk(integrand, -Inf, Inf; rtol = 1e-8, atol = 1e-10)
# 	return result * 2
# end
# function f(q, α)
# 	return 2π * erfcx(q / (2 * α))/q
# end
