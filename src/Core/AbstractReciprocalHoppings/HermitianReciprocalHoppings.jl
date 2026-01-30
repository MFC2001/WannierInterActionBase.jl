"""
	HermitianReciprocalHoppings{T <: Number, U <: AbstractReciprocalHoppings{T}} <: AbstractReciprocalHoppings{T}

```julia
julia> hr isa wannier90_hr
true

julia> hrh = HermitianReciprocalHoppings(hr);

julia> k isa ReducedCoordinates
true

julia> hrh(k) isa Hermitian
true

```
"""
struct HermitianReciprocalHoppings{T <: Number, U <: AbstractReciprocalHoppings{T}} <: AbstractReciprocalHoppings{T}
	core::U
	uplo::Symbol
end

@inline function (hrh::HermitianReciprocalHoppings)(k::ReducedCoordinates)
	A = Matrix{ComplexF64}(undef, numorb(hrh), numorb(hrh))
	return hrh(A, k)
end
@inline function (hrh::HermitianReciprocalHoppings)(A::AbstractMatrix, k::ReducedCoordinates)
	size(A) == (numorb(hrh), numorb(hrh)) || error("Buffer size mismatch.")
	_hrh_triangular_term!(A, Val(hrh.uplo), hrh, k)
	return Hermitian(A, hrh.uplo)
end

@inline function (hrh::HermitianReciprocalHoppings)(k::ReducedCoordinates, orblocat)
	A = Matrix{ComplexF64}(undef, numorb(hrh), numorb(hrh))
	return hrh(A, k, orblocat)
end
@inline function (hrh::HermitianReciprocalHoppings)(A::AbstractMatrix, k::ReducedCoordinates, orblocat)
	size(A) == (numorb(hrh), numorb(hrh)) || error("Buffer size mismatch.")
	length(orblocat) == numorb(hrh) || error("Mismatched orbital locations.")
	_hrh_triangular_term!(A, Val(hrh.uplo), hrh, k, orblocat)
	return Hermitian(A, hrh.uplo)
end

@inline function _hrh_triangular_term!(A, ::Val{:U}, hrh, parameters...)
	@inbounds for j in 1:numorb(hrh), i in 1:j
		calculate_value = _Arh_getindex(hrh, i, j, parameters...)
		A[i, j] = calculate_value
		A[j, i] = conj(calculate_value)
	end
end
@inline function _hrh_triangular_term!(A, ::Val{:L}, hrh, parameters...)
	@inbounds for j in 1:numorb(hrh), i in j:numorb(hrh)
		calculate_value = _Arh_getindex(hrh, i, j, parameters...)
		A[i, j] = calculate_value
		A[j, i] = conj(calculate_value)
	end
end


@inline function (hrh::HermitianReciprocalHoppings)(::Val{:partial}, lattice::Lattice, k::ReducedCoordinates, orblocat)
	A = Array{ComplexF64}(undef, numorb(hrh), numorb(hrh), 3)
	return hrh(A, Val(:partial), lattice, k, orblocat)
end
@inline function (hrh::HermitianReciprocalHoppings)(::Val{:partial}, A::AbstractArray, lattice::Lattice, k::ReducedCoordinates, orblocat)
	size(A) == (numorb(hrh), numorb(hrh), 3) || error("Buffer size mismatch.")
	length(orblocat) == numorb(hrh) || error("Mismatched orbital locations.")
	_hrh_triangular_term!(A, Val(hrh.uplo), hrh, Val(:partial), lattice, k, orblocat)
	return A
end

@inline function _hrh_triangular_term!(A, ::Val{:U}, hrh, ::Val{:partial}, lattice, k, orblocat)
	@inbounds for j in 1:numorb(hrh), i in 1:j
		calculate_value = _Arh_getindex(hrh, i, j, Val(:partial), lattice, k, orblocat)
		A[i, j, :] .= calculate_value
		A[j, i, :] .= conj.(calculate_value)
	end
end
@inline function _hrh_triangular_term!(A, ::Val{:L}, hrh, ::Val{:partial}, lattice, k, orblocat)
	@inbounds for j in 1:numorb(hrh), i in j:numorb(hrh)
		calculate_value = _Arh_getindex(hrh, i, j, Val(:partial), lattice, k, orblocat)
		A[i, j, :] .= calculate_value
		A[j, i, :] .= conj.(calculate_value)
	end
end



function Base.show(io::IO, hrh::HermitianReciprocalHoppings)
	print(io, "HermitianReciprocalHoppings with $(sum(Nhop(hrh.core))) hoppings and $(numorb(hrh.core)) orbitals.")
end

function HermitianReciprocalHoppings(hr::wannier90_hr{T}, uplo::Symbol = :U) where {T}
	if uplo ∉ (:U, :L)
		throw(ArgumentError("uplo must be either :U or :L, got :$uplo"))
	end
	rh = BaseReciprocalHoppings(hr)
	return HermitianReciprocalHoppings(rh, uplo)
end
function HermitianReciprocalHoppings(rh::BaseReciprocalHoppings{T}, uplo::Symbol = :U) where {T}
	if uplo ∉ (:U, :L)
		throw(ArgumentError("uplo must be either :U or :L, got :$uplo"))
	end
	return HermitianReciprocalHoppings{T, BaseReciprocalHoppings{T}}(rh, uplo)
end
function HermitianReciprocalHoppings(rh::ReciprocalHoppings{T, U}, uplo::Symbol = :U) where {T, U}
	if uplo ∉ (:U, :L)
		throw(ArgumentError("uplo must be either :U or :L, got :$uplo"))
	end
	return HermitianReciprocalHoppings{T, U}(rh.core, uplo)
end
function HermitianReciprocalHoppings(rh::HermitianReciprocalHoppings{T, U}, uplo::Symbol = :U) where {T, U}
	if uplo ∉ (:U, :L)
		throw(ArgumentError("uplo must be either :U or :L, got :$uplo"))
	end
	return HermitianReciprocalHoppings{T, U}(rh.core, uplo)
end

Base.zero(::HermitianReciprocalHoppings{T, U}, norb::Integer) where {T, U} =
	HermitianReciprocalHoppings{T, U}(zero(U, norb))
Base.zero(::Type{HermitianReciprocalHoppings{T, U}}, norb::Integer) where {T, U} =
	HermitianReciprocalHoppings{T, U}(zero(U, norb))

function Base.union(rh₁::HermitianReciprocalHoppings{T₁, U₁}, rh₂::AbstractReciprocalHoppings{T₂}) where {T₁, U₁, T₂}
	return HermitianReciprocalHoppings(union(rh₁.core, BaseReciprocalHoppings(rh₂)), rh₁.uplo)
end

function spinrh(hrh::HermitianReciprocalHoppings{T, U}; mode = conj) where {T, U}
	return HermitianReciprocalHoppings(spinrh(hrh.core; mode), hrh.uplo)
end
