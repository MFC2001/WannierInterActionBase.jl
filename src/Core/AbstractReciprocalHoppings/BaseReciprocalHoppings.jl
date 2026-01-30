struct BaseReciprocalHoppings{T} <: AbstractReciprocalHoppings{T}
	norb::Int
	Nhop::Matrix{Int}
	hops::Matrix{Vector{Hopping{T}}}
end
@inline numorb(brh::BaseReciprocalHoppings) = brh.norb
@inline Nhop(brh::BaseReciprocalHoppings) = brh.Nhop
@inline hops(brh::BaseReciprocalHoppings) = brh.hops
@inline BaseReciprocalHoppings(brh::BaseReciprocalHoppings) = brh

@inline function _Arh_getindex(brh::BaseReciprocalHoppings, I::CartesianIndex{2}, k::ReducedCoordinates)
	@inbounds iszero(brh.Nhop[I]) ? 0 :
			  sum(hop -> hop.t * hopphase(hop, k), brh.hops[I])
end
@inline function _Arh_getindex(brh::BaseReciprocalHoppings, I::CartesianIndex{2}, k::ReducedCoordinates, orblocat)
	(i, j) = Tuple(I)
	@inbounds iszero(brh.Nhop[i, j]) ? 0 :
			  cis(2π * (k ⋅ (orblocat[j] - orblocat[i]))) * sum(hop -> hop.t * hopphase(hop, k), brh.hops[i, j])
end
@inline function _Arh_getindex(brh::BaseReciprocalHoppings, I::CartesianIndex{2}, ::Val{:partial}, lattice::Lattice, k::ReducedCoordinates, orblocat)
	(i, j) = Tuple(I)
	@inbounds dorb = orblocat[j] - orblocat[i]
	@inbounds iszero(brh.Nhop[i, j]) ? zeros(ComplexF64, 3) :
			  lattice * (im * cis(2π * (k ⋅ dorb)) * sum(hop -> hop.t * hopphase(hop, k) * (hop.R + dorb), brh.hops[i, j]))
end
@inline function _Arh_getindex(brh::BaseReciprocalHoppings, i::Integer, j::Integer, k::ReducedCoordinates)
	@inbounds iszero(brh.Nhop[i, j]) ? 0 :
			  sum(hop -> hop.t * hopphase(hop, k), brh.hops[i, j])
end
@inline function _Arh_getindex(brh::BaseReciprocalHoppings, i::Integer, j::Integer, k::ReducedCoordinates, orblocat)
	@inbounds iszero(brh.Nhop[i, j]) ? 0 :
			  cis(2π * (k ⋅ (orblocat[j] - orblocat[i]))) * sum(hop -> hop.t * hopphase(hop, k), brh.hops[i, j])
end
@inline function _Arh_getindex(brh::BaseReciprocalHoppings, i::Integer, j::Integer, ::Val{:partial}, lattice::Lattice, k::ReducedCoordinates, orblocat)
	@inbounds dorb = orblocat[j] - orblocat[i]
	@inbounds iszero(brh.Nhop[i, j]) ? zeros(ComplexF64, 3) :
			  lattice * (im * cis(2π * (k ⋅ dorb)) * sum(hop -> hop.t * hopphase(hop, k) * (hop.R + dorb), brh.hops[i, j]))
end


function BaseReciprocalHoppings(hr::wannier90_hr{T}) where {T}
	norb = numorb(hr)
	if sort(hr.orbindex) == 1:norb
		hops = _buildhops_from_hr(norb, hr.path, hr.value)
		Nhop = length.(hops)
	else
		error("Won't construct `ReciprocalHoppings from a nonstandard `wannier90_hr`.")
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

Base.zero(::BaseReciprocalHoppings{T}, norb::Integer) where {T <: Number} = zero(BaseReciprocalHoppings{T}, norb)
function Base.zero(::Type{BaseReciprocalHoppings{T}}, norb::Integer) where {T <: Number}
	hops = [Vector{Hopping{T}}(undef, 0) for _ in CartesianIndices((norb, norb))]
	Nhop = zeros(Int, norb, norb)
	return BaseReciprocalHoppings{T}(norb, hops, Nhop)
end
function Base.convert(::BaseReciprocalHoppings{T₁}, brh::BaseReciprocalHoppings{T₂}) where {T₁, T₂}
	return BaseReciprocalHoppings{T₁}(brh.norb, brh.Nhop, brh.hops)
end
function Base.convert(::Type{BaseReciprocalHoppings{T₁}}, brh::BaseReciprocalHoppings{T₂}) where {T₁, T₂}
	return BaseReciprocalHoppings{T₁}(brh.norb, brh.Nhop, brh.hops)
end
function Base.union(brh₁::BaseReciprocalHoppings{T₁}, brh₂::BaseReciprocalHoppings{T₂}) where {T₁, T₂}
	numorb(brh₁) == numorb(brh₂) || error("Mismatched BaseReciprocalHoppings.")

	norb = numorb(brh₁)
	T = promote_type(T₁, T₂)
	hops_M = [Vector{Hopping{T}}(undef, 0) for _ in CartesianIndices((norb, norb))]
	for I in CartesianIndices((norb, norb))
		hops_M[I] = _unique_hops_ij_rh!([hops(brh₁)[I]; hops(brh₂)[I]])
	end
	Nhop_M = length.(hops_M)

	return BaseReciprocalHoppings{T}(norb, Nhop_M, hops_M)
end
function _union_hops_ij_rh!(hops₁::AbstractVector{Hopping{T₁}}, hops₂::AbstractVector{Hopping{T₂}}) where {T₁, T₂}
	append!(hops₁, hops₂)
	return _unique_hops_ij_rh!(hops₁)
end
function _unique_hops_ij_rh!(hops::AbstractVector{Hopping{T}}) where {T}
	if isempty(hops)
		return hops
	end
	unique_R = Vector{ReducedCoordinates{Int}}(undef, 0)
	n = 1
	while true
		R = path(hops[n])
		i = findfirst(isequal(R), unique_R)
		if isnothing(i)
			n += 1
			push!(unique_R, R)
		else
			hops[i] = similar(hops[i], hops[i].t + hops[n].t)
			deleteat!(hops, n)
		end
		if n > length(hops)
			break
		end
	end
	return hops
end

function wannier90_hr(brh::BaseReciprocalHoppings{T}) where {T}
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

	return wannier90_hr(path, value; orbindex = collect(1:numorb(brh)), hrsort = true)
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
