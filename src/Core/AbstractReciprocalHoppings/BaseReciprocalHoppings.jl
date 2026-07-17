struct BaseReciprocalHoppings{T} <: AbstractReciprocalHoppings{T}
	norb::Int
	Nhop::Matrix{Int}
	hops::Matrix{Vector{Hopping{T}}}
end
@inline numorb(brh::BaseReciprocalHoppings) = brh.norb
@inline Nhop(brh::BaseReciprocalHoppings) = brh.Nhop
@inline hops(brh::BaseReciprocalHoppings) = brh.hops
@inline BaseReciprocalHoppings(brh::BaseReciprocalHoppings) = brh

# @inline function _Arh_getindex(brh::BaseReciprocalHoppings, I::CartesianIndex{2}, k::ReducedCoordinates)
# 	@inbounds iszero(brh.Nhop[I]) ? 0 :
# 			  sum(hop -> hop.t * hopphase(hop, k), brh.hops[I])
# end
# @inline function _Arh_getindex(brh::BaseReciprocalHoppings, I::CartesianIndex{2}, k::ReducedCoordinates, orblocat::AbstractVector{<:ReducedCoordinates})
# 	(i, j) = Tuple(I)
# 	@inbounds iszero(brh.Nhop[i, j]) ? 0 :
# 			  cispi(2 * (k ⋅ (orblocat[j] - orblocat[i]))) * sum(hop -> hop.t * hopphase(hop, k), brh.hops[i, j])
# end
# @inline function _Arh_getindex(brh::BaseReciprocalHoppings, I::CartesianIndex{2}, ::Val{:partial}, lattice::Lattice, k::ReducedCoordinates, orblocat::AbstractVector{<:ReducedCoordinates})
# 	(i, j) = Tuple(I)
# 	@inbounds dorb = orblocat[j] - orblocat[i]
# 	@inbounds iszero(brh.Nhop[i, j]) ? zeros(ComplexF64, 3) :
# 			  lattice * (im * cispi(2 * (k ⋅ dorb)) * sum(hop -> hop.t * hopphase(hop, k) * (hop.R + dorb), brh.hops[i, j]))
# end
# @inline function _Arh_getindex(brh::BaseReciprocalHoppings, i::Integer, j::Integer, k::ReducedCoordinates)
# 	@inbounds iszero(brh.Nhop[i, j]) ? 0 :
# 			  sum(hop -> hop.t * hopphase(hop, k), brh.hops[i, j])
# end
# @inline function _Arh_getindex(brh::BaseReciprocalHoppings, i::Integer, j::Integer, k::ReducedCoordinates, orblocat::AbstractVector{<:ReducedCoordinates})
# 	@inbounds iszero(brh.Nhop[i, j]) ? 0 :
# 			  cispi(2 * (k ⋅ (orblocat[j] - orblocat[i]))) * sum(hop -> hop.t * hopphase(hop, k), brh.hops[i, j])
# end
# @inline function _Arh_getindex(brh::BaseReciprocalHoppings, i::Integer, j::Integer, ::Val{:partial}, lattice::Lattice, k::ReducedCoordinates, orblocat::AbstractVector{<:ReducedCoordinates})
# 	@inbounds dorb = orblocat[j] - orblocat[i]
# 	@inbounds iszero(brh.Nhop[i, j]) ? zeros(ComplexF64, 3) :
# 			  lattice * (im * cispi(2 * (k ⋅ dorb)) * sum(hop -> hop.t * hopphase(hop, k) * (hop.R + dorb), brh.hops[i, j]))
# end

@inline function _Arh_calc!(A::AbstractMatrix, brh::BaseReciprocalHoppings, ij::Vector{Tuple{Int, Int}}, k::ReducedCoordinates)
	@inbounds for (i, j) in ij
		A[i, j] = iszero(brh.Nhop[i, j]) ? 0 :
				  sum(hop -> hop.t * hopphase(hop, k), brh.hops[i, j])
	end
	return A
end
@inline function _Arh_calc!(A::AbstractMatrix, brh::BaseReciprocalHoppings, ij::Vector{Tuple{Int, Int}}, k::ReducedCoordinates, orblocat::AbstractVector{<:ReducedCoordinates})
	@inbounds for (i, j) in ij
		A[i, j] = iszero(brh.Nhop[i, j]) ? 0 :
				  cispi(2 * (k ⋅ (orblocat[j] - orblocat[i]))) * sum(hop -> hop.t * hopphase(hop, k), brh.hops[i, j])
	end
	return A
end
@inline function _Arh_calc!(A::AbstractArray, brh::BaseReciprocalHoppings, ij::Vector{Tuple{Int, Int}}, ::Val{:partial}, k::ReducedCoordinates, lattice::Lattice, orblocat::AbstractVector{<:ReducedCoordinates})
	@inbounds for (i, j) in ij
		dorb = orblocat[j] - orblocat[i]
		if iszero(brh.Nhop[i, j])
			A[i, j, :] .= zero(ComplexF64)
		else
			A[i, j, :] .= lattice * (im * cispi(2 * (k ⋅ dorb)) * sum(hop -> hop.t * hopphase(hop, k) * (hop.R + dorb), brh.hops[i, j]))
		end
	end
	return A
end

function Base.filter(F, brh::BaseReciprocalHoppings{T}) where {T}
	norb = numorb(brh)
	hops = Matrix{Vector{Hopping{T}}}(undef, norb, norb)
	for j in 1:norb, i in 1:norb
		hops[i, j] = filter(F, brh.hops[i, j])
	end
	Nhop = length.(hops)
	return BaseReciprocalHoppings{T}(norb, Nhop, hops)
end
function Base.filter!(F, brh::BaseReciprocalHoppings{T}) where {T}
	norb = numorb(brh)
	for j in 1:norb, i in 1:norb
		filter!(F, brh.hops[i, j])
	end
	brh.Nhop .= length.(brh.hops)
	return brh
end
function prune(brh::BaseReciprocalHoppings{T}, cutoff::Real) where {T <: Real}
	@assert cutoff ≥ 0 "cutoff should be positive!"
	f(hop) = abs(value(hop)) ≥ cutoff
	return filter(f, brh)
end
function prune(brh::BaseReciprocalHoppings{T}, cutoff::Real) where {T <: Complex}
	@assert cutoff ≥ 0 "cutoff should be positive!"
	cutoff2 = cutoff^2
	f(hop) = abs2(value(hop)) ≥ cutoff2
	return filter(f, brh)
end
function prune!(brh::BaseReciprocalHoppings{T}, cutoff::Real) where {T <: Real}
	@assert cutoff ≥ 0 "cutoff should be positive!"
	f(hop) = abs(value(hop)) ≥ cutoff
	return filter!(f, brh)
end
function prune!(brh::BaseReciprocalHoppings{T}, cutoff::Real) where {T <: Complex}
	@assert cutoff ≥ 0 "cutoff should be positive!"
	cutoff2 = cutoff^2
	f(hop) = abs2(value(hop)) ≥ cutoff2
	return filter!(f, brh)
end
function Base.real(brh::BaseReciprocalHoppings{T}) where {T <: Complex}
	norb = numorb(brh)
	realT = real(T)
	hops = Matrix{Vector{Hopping{realT}}}(undef, norb, norb)
	for j in 1:norb, i in 1:norb
		hops[i, j] = real.(brh.hops[i, j])
	end
	return BaseReciprocalHoppings{realT}(norb, copy(brh.Nhop), hops)
end
Base.real(::Type{BaseReciprocalHoppings{T}}) where {T <: Real} = BaseReciprocalHoppings{T}
Base.real(::Type{BaseReciprocalHoppings{T}}) where {T <: Complex} = BaseReciprocalHoppings{real(T)}
function Base.isreal(brh::BaseReciprocalHoppings{T}) where {T <: Complex}
	result = true
	norb = numorb(brh)
	for j in 1:norb, i in 1:norb
		result &= all(isreal, brh.hops[i, j])
	end
	return result
end
function Base.iszero(brh::BaseReciprocalHoppings)
	result = true
	norb = numorb(brh)
	for j in 1:norb, i in 1:norb
		result &= (iszero(brh.Nhop[i, j]) || all(iszero, brh.hops[i, j]))
	end
	return result
end
@inline Base.zero(::BaseReciprocalHoppings{T}, norb::Integer) where {T} = zero(BaseReciprocalHoppings{T}, norb)
function Base.zero(::Type{BaseReciprocalHoppings{T}}, norb::Integer) where {T <: Number}
	hops = [Vector{Hopping{T}}(undef, 0) for _ in CartesianIndices((norb, norb))]
	Nhop = zeros(Int, norb, norb)
	return BaseReciprocalHoppings{T}(norb, hops, Nhop)
end
function Base.convert(::Type{BaseReciprocalHoppings{T₁}}, brh::BaseReciprocalHoppings{T₂}) where {T₁, T₂}
	return BaseReciprocalHoppings{T₁}(brh.norb, brh.Nhop, brh.hops)
end

"""
	translate(brh::BaseReciprocalHoppings{T}, centres::Pair{<:Integer, <:AbstractVector{<:Integer}}...)

You only need to input the wannier centres R which isn't in the [0, 0, 0] cell.

See also `translate!`.
"""
function translate(brh::BaseReciprocalHoppings{T}, aim_centres::Pair{<:Integer, <:AbstractVector{<:Integer}}...) where {T}
	norb = numorb(brh)
	centres = [ReducedCoordinates{Int}(0, 0, 0) for _ in Base.OneTo(norb)]
	for (iw, R) in aim_centres
		centres[iw] = ReducedCoordinates{Int}(R)
	end
	new_hops = Matrix{Vector{Hopping{T}}}(undef, norb, norb)
	for j in Base.OneTo(norb), i in Base.OneTo(norb)
		Nhop_ij = Nhop(brh)[i, j]
		if iszero(Nhop_ij)
			new_hops[i, j] = Vector{Hopping{T}}(undef, 0)
			continue
		end
		dR = centres[j] - centres[i]
		hops_ij = hops(brh)[i, j]
		if iszero(dR)
			new_hops[i, j] = copy(hops_ij)
		else
			new_hops_ij = Vector{Hopping{T}}(undef, Nhop_ij)
			for ihop in Base.OneTo(Nhop_ij)
				hop = hops_ij[ihop]
				new_hops_ij[ihop] = Hopping{T}(i, j, hop.R + dR, hop.t)
			end
		end
	end
	return BaseReciprocalHoppings{T}(norb, new_hops, copy(Nhop(brh)))
end

@inline _brh_type(::BaseReciprocalHoppings{T}) where {T} = T
function Base.merge(brh₁::BaseReciprocalHoppings{T₁}, brh₂::BaseReciprocalHoppings{T₂}) where {T₁, T₂}
	# 不能定义成brhs...形式。
	norb = numorb(brh₁)
	norb == numorb(brh₂) || error("Mismatched BaseReciprocalHoppings.")

	T = promote_type(T₁, T₂)

	all_hops = Matrix{Vector{Hopping{T}}}(undef, norb, norb)
	for I in CartesianIndices((norb, norb))
		all_hops[I] = merge_ijhops(hops(brh₁)[I], hops(brh₂)[I])
	end
	all_Nhop = length.(all_hops)

	return BaseReciprocalHoppings{T}(norb, all_Nhop, all_hops)
end
function Base.merge(brh₁::BaseReciprocalHoppings{T₁}, hr₂::HR{T₂}) where {T₁, T₂}
	norb = numorb(brh₁)
	norb == numorb(hr₂) || error("Mismatched BaseReciprocalHoppings.")

	return merge(brh₁, BaseReciprocalHoppings(hr₂))
end
function Base.merge(brh::BaseReciprocalHoppings{T}, arrays::Union{AbstractArray{<:Hopping}, Hopping}...) where {T}

	newT = reduce(promote_type, _hop_amplitude_type(arr) for arr in arrays; init = T)

	norb = numorb(brh)
	hops_new = [Vector{Hopping{newT}}(undef, 0) for _ in CartesianIndices((norb, norb))]
	arrays = Iterators.flatten(arrays)
	for hop in arrays
		push!(hops_new[hop.i, hop.j], hop)
	end

	all_hops = Matrix{Vector{Hopping{newT}}}(undef, norb, norb)
	for I in CartesianIndices((norb, norb))
		all_hops[I] = merge_ijhops(hops(brh)[I], hops_new[I])
	end
	all_Nhop = length.(all_hops)

	return BaseReciprocalHoppings{newT}(norb, all_Nhop, all_hops)
end
function Base.merge!(brh₁::BaseReciprocalHoppings{T₁}, brh₂::BaseReciprocalHoppings{T₂}) where {T₁, T₂}
	# 不能定义成brhs...形式。
	norb = numorb(brh₁)
	norb == numorb(brh₂) || error("Mismatched BaseReciprocalHoppings.")

	promote_type(T₁, T₂) == T₁ || error("Incompatible types for `merge!`")

	for I in CartesianIndices((norb, norb))
		merge_ijhops!(hops(brh₁)[I], hops(brh₂)[I])
	end
	Nhop(brh₁) .= length.(hops(brh₁))

	return brh₁
end
function Base.merge!(brh₁::BaseReciprocalHoppings{T₁}, hr₂::HR{T₂}) where {T₁, T₂}
	norb = numorb(brh₁)
	norb == numorb(hr₂) || error("Mismatched BaseReciprocalHoppings.")

	return merge!(brh₁, BaseReciprocalHoppings(hr₂))
end
function Base.merge!(brh::BaseReciprocalHoppings{T}, arrays::Union{AbstractArray{<:Hopping}, Hopping}...) where {T}

	reduce(promote_type, _hop_amplitude_type(arr) for arr in arrays; init = T) == T || error("Incompatible types for `merge!`")

	norb = numorb(brh)
	hops_new = [Vector{Hopping{T}}(undef, 0) for _ in CartesianIndices((norb, norb))]
	arrays = Iterators.flatten(arrays)
	for hop in arrays
		push!(hops_new[hop.i, hop.j], hop)
	end

	for I in CartesianIndices((norb, norb))
		merge_ijhops!(hops(brh)[I], hops_new[I])
	end
	Nhop(brh) .= length.(hops(brh))

	return brh
end

function spinrh(brh::BaseReciprocalHoppings{T}; mode = conj) where {T}
	norb = numorb(brh)
	hops_up = deepcopy(hops(brh))
	hops_dn = map(hops_up) do hops
		map(hop -> Hopping{T}(hop.i, hop.j, hop.R, mode(hop.t)), hops)
	end
	hops_0 = [Vector{Hopping{T}}(undef, 0) for _ in CartesianIndices((norb, norb))]
	spinhops = [
		hops_up hops_0;
		hops_0 hops_dn
	]
	spinNhop = length.(spinhops)
	return BaseReciprocalHoppings(norb * 2, spinNhop, spinhops)
end

function BaseReciprocalHoppings(hr::HR{T}) where {T}
	norb = numorb(hr)
	if sort(hr.orbindex) == 1:norb
		hops = _buildhops_from_hr(norb, hr.path, hr.value)
		Nhop = length.(hops)
	else
		error("Won't construct `ReciprocalHoppings from a nonstandard `HR`.")
	end
	return BaseReciprocalHoppings{T}(norb, Nhop, hops)
end
function _buildhops_from_hr(norb, path::AbstractMatrix, value::AbstractVector{T}) where {T}
	hops = [Vector{Hopping{T}}(undef, 0) for _ in CartesianIndices((norb, norb))]
	for i in axes(path, 2)
		push!(hops[path[4, i], path[5, i]], Hopping(path[:, i], value[i]))
	end
	return hops
end
function BaseReciprocalHoppings(int::WanIntijklRFull{T}) where {T}
	np = length(int.pair)
	hops = [Vector{Hopping{T}}(undef, 0) for _ in CartesianIndices((np, np))]
	nR = length(int.R)
	for ip2 in Base.OneTo(np), ip1 in Base.OneTo(np)
		sizehint!(hops[ip1, ip2], nR)
	end
	for (ip1, ip2, iR, v) in enumerate(int)
		R = int.R[iR]
		push!(hops[ip1, ip2], Hopping{T}(ip1, ip2, R, v))
	end
	Nhop = length.(hops)
	return BaseReciprocalHoppings{T}(np, Nhop, hops)
end
function BaseReciprocalHoppings(int::WanIntijklRHermitian{T}) where {T}
	np = length(int.pair)
	hops = [Vector{Hopping{T}}(undef, 0) for _ in CartesianIndices((np, np))]
	nR = length(int.R)
	for ip2 in Base.OneTo(np), ip1 in Base.OneTo(np)
		sizehint!(hops[ip1, ip2], nR)
	end
	for (ip1, ip2, iR, iconjR, v) in zip(int.ips, int.jps, int.iRs, int.iconjRs, int.vals)
		R = int.R[iR]
		push!(hops[ip1, ip2], Hopping{T}(ip1, ip2, R, v))
		R = int.R[iconjR]
		ip1 = int.conjpair[ip1]
		ip2 = int.conjpair[ip2]
		push!(hops[ip1, ip2], Hopping{T}(ip1, ip2, R, conj(v)))
	end
	Nhop = length.(hops)
	return BaseReciprocalHoppings{T}(np, Nhop, hops)
end

function HR(brh::BaseReciprocalHoppings{T}) where {T}
	N = sum(Nhop(brh))
	path = Matrix{Int}(undef, 5, N)
	value = Vector{T}(undef, N)

	n = 0
	for I in CartesianIndices((numorb(brh), numorb(brh)))
		for hop in hops(brh)[I]
			n += 1
			path[:, n] = [hop.R; hop.i; hop.j]
			value[n] = hop.t
		end
	end

	return HR(path, value; orbindex = collect(1:numorb(brh)), hrsort = true)
end
