
struct ΓReciprocalHoppings{T} <: ReciprocalHoppings{T} 
	core::U
end

@inline function (rh::ΓReciprocalHoppings{T})(paras...) where {T}
	A = Matrix{T}(undef, rh.norb, rh.norb)
	return rh(A, paras...)
end
# @inline function (rh::ΓReciprocalHoppings{T})(k::AbstractVector) where {T}
# 	A = Matrix{T}(undef, rh.norb, rh.norb)
# 	return rh(A, k)
# end
# @inline function (rh::ΓReciprocalHoppings{T})(k::AbstractVector, orblocat) where {T}
# 	A = Matrix{T}(undef, rh.norb, rh.norb)
# 	return rh(A, k, orblocat)
# end
@inline function (rh::ΓReciprocalHoppings{T})(I::CartesianIndex{2}, paras...) where {T}
	iszero(rh.Nhop[I]) ? 0 :
	sum(hop -> hop.t, rh.hops[I])
end
# @inline function (rh::ΓReciprocalHoppings{T})(I::CartesianIndex{2}, k) where {T}
# 	iszero(rh.Nhop[I]) ? 0 :
# 	sum(hop -> hop.t, rh.hops[I])
# end
# @inline function (rh::ΓReciprocalHoppings{T})(I::CartesianIndex{2}, k, orblocat) where {T}
# 	(i, j) = Tuple(I)
# 	iszero(rh.Nhop[i, j]) ? 0 :
# 	sum(hop -> hop.t, rh.hops[i, j])
# end
@inline function (rh::ΓReciprocalHoppings{T})(i::Integer, j::Integer, paras...) where {T}
	iszero(rh.Nhop[i, j]) ? 0 :
	sum(hop -> hop.t, rh.hops[i, j])
end
# @inline function (rh::ΓReciprocalHoppings{T})(i::Integer, j::Integer, k) where {T}
# 	iszero(rh.Nhop[i, j]) ? 0 :
# 	sum(hop -> hop.t, rh.hops[i, j])
# end
# @inline function (rh::ΓReciprocalHoppings{T})(i::Integer, j::Integer, k, orblocat) where {T}
# 	iszero(rh.Nhop[i, j]) ? 0 :
# 	sum(hop -> hop.t, rh.hops[i, j])
# end

#这里可以将其分为两种，实数版和复数版
