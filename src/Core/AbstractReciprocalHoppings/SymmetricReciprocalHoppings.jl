abstract type HopMap end

struct SymmetricReciprocalHoppings{T} <: ReciprocalHoppings{T}
	norb::Int
	Nhop::Matrix{Int}
	hops::Matrix{Vector{Hopping{T}}}
	irredhop::Vector{Hopping{T}}
	hopmap::Matrix{Vector{HopMap}}
end
#为每种映射定义一个数据类型，它们均为抽象类型HopMap的子类型，并且均支持行为hopmap(irredhop).
#其内部可能会包含一系列数据，涉及到计算系数等。

function refresh!(srh::SymmetricReciprocalHoppings)
	for I in CartesianIndices(srh.hops)
		srh.hops[I] = map(f -> f(srh.value), srh.hopmap[I])
	end
end