
function _RealCoulombKeldysh(::Val{:Struve}, v, r)::Float64
	r = r / v.r₀
	return v.coff * (Struve.struveh(0, r) - bessely0(r)) * π / 2
end
