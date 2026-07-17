export WilsonLoop, WilsonLoopWaveGetter, Index_VecOrMat
abstract type WilsonLoopWaveGetter{Nstates} end
# function AbstractWaveGetter(kpoint) end

# function WilsonLoop(vertex::AbstractVector)
	
# end

# 这里只是计算了一个路径，如果需要路径是环路，需要loop[end] == loop[1]
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

struct IndexToU{N, T <: AbstractArray} <: WilsonLoopWaveGetter{N}
	U::Vector{T}
end
(u::IndexToU)(k) = u.U[k]

function WilsonLoopWaveGetter(wave::AbstractVector{T}) where {T <: AbstractArray}
	nband = size.(wave, 2)
	if any(x -> x ≠ nband[1], nband)
		error("Wrong wave!")
	end
	return IndexToU{nband[1], T}(collect(wave))
end
