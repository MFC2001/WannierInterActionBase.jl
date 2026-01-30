abstract type WilsonLoopWaveGetter{Nstates} end
# function AbstractWaveGetter(kpoint) end

function WilsonLoop(loop::AbstractVector, u_getter::WilsonLoopWaveGetter{1})
	u1 = u_getter(loop[1])
	W = 1
	for k in loop[2:end]
		u2 = u_getter(k)
		uu = u2 ⋅ u1
		W *= uu
		u1 = u2
	end
	return W / abs(W)
end
function WilsonLoop(loop::AbstractVector, u_getter::WilsonLoopWaveGetter)
	u1 = u_getter(loop[1])
	W = 1
	for k in loop[2:end]
		u2 = u_getter(k)
		uu = u2' * u1
		F = svd(uu)
		uu = F.U * F.Vt
		W = uu * W
		u1 = u2
	end
	return W
end

struct Index_VecOrMat{N} <: WilsonLoopWaveGetter{N}
	U::Vector{VecOrMat{ComplexF64}}
end
(u::Index_VecOrMat)(k) = u.U[k]
