"""
	HermitianReciprocalHoppings{T <: Number, U <: AbstractReciprocalHoppings{T}} <: AbstractReciprocalHoppings{T}

```julia
julia> hr isa AbstractHR
true

julia> hrh = HermitianReciprocalHoppings(hr);

julia> k isa ReducedCoordinates
true

julia> hrh(k) isa Hermitian
true

```
"""
struct HermitianReciprocalHoppings{T<:Number,U<:AbstractReciprocalHoppings{T}} <: AbstractReciprocalHoppings{T}
    core::U
    uplo::Symbol
    ij::Vector{Tuple{Int,Int}}
end

@inline function (hrh::HermitianReciprocalHoppings)(k::ReducedCoordinates)
    A = Matrix{ComplexF64}(undef, numorb(hrh), numorb(hrh))
    _Arh_calc!(A, hrh.core, hrh.ij, k)
    _hrh_fill!(A, numorb(hrh), Val(hrh.uplo))
    return Hermitian(A, hrh.uplo)
end
@inline function (hrh::HermitianReciprocalHoppings)(A::AbstractMatrix, k::ReducedCoordinates)
    size(A) == (numorb(hrh), numorb(hrh)) || error("Buffer size mismatch.")
    _Arh_calc!(A, hrh.core, hrh.ij, k)
    _hrh_fill!(A, numorb(hrh), Val(hrh.uplo))
    return Hermitian(A, hrh.uplo)
end

@inline function (hrh::HermitianReciprocalHoppings)(::Val{:partial}, k::ReducedCoordinates, lattice::Lattice)
    A = Array{ComplexF64}(undef, numorb(hrh), numorb(hrh), 3)
    _Arh_calc!(A, hrh.core, hrh.ij, Val(:partial), k, lattice)
    for i in 1:3
        _hrh_fill!(view(A, :, :, i), numorb(hrh), Val(hrh.uplo))
    end
    return A
end
@inline function (hrh::HermitianReciprocalHoppings)(::Val{:partial}, A::AbstractArray, k::ReducedCoordinates, lattice::Lattice)
    size(A) == (numorb(hrh), numorb(hrh), 3) || error("Buffer size mismatch.")
    _Arh_calc!(A, hrh.core, hrh.ij, Val(:partial), k, lattice)
    for i in 1:3
        _hrh_fill!(view(A, :, :, i), numorb(hrh), Val(hrh.uplo))
    end
    return A
end

@inline function (hrh::HermitianReciprocalHoppings)(k::ReducedCoordinates, orblocat::AbstractVector{<:ReducedCoordinates})
    A = Matrix{ComplexF64}(undef, numorb(hrh), numorb(hrh))
    length(orblocat) == numorb(hrh) || error("Mismatched orbital locations.")
    _Arh_calc!(A, hrh.core, hrh.ij, k, orblocat)
    _hrh_fill!(A, numorb(hrh), Val(hrh.uplo))
    return Hermitian(A, hrh.uplo)
end
@inline function (hrh::HermitianReciprocalHoppings)(A::AbstractMatrix, k::ReducedCoordinates, orblocat::AbstractVector{<:ReducedCoordinates})
    size(A) == (numorb(hrh), numorb(hrh)) || error("Buffer size mismatch.")
    length(orblocat) == numorb(hrh) || error("Mismatched orbital locations.")
    _Arh_calc!(A, hrh.core, hrh.ij, k, orblocat)
    _hrh_fill!(A, numorb(hrh), Val(hrh.uplo))
    return Hermitian(A, hrh.uplo)
end

@inline function (hrh::HermitianReciprocalHoppings)(::Val{:partial}, k::ReducedCoordinates, lattice::Lattice, orblocat::AbstractVector{<:ReducedCoordinates})
    A = Array{ComplexF64}(undef, numorb(hrh), numorb(hrh), 3)
    length(orblocat) == numorb(hrh) || error("Mismatched orbital locations.")
    _Arh_calc!(A, hrh.core, hrh.ij, Val(:partial), k, lattice, orblocat)
    for i in 1:3
        _hrh_fill!(view(A, :, :, i), numorb(hrh), Val(hrh.uplo))
    end
    return A
end
@inline function (hrh::HermitianReciprocalHoppings)(::Val{:partial}, A::AbstractArray, k::ReducedCoordinates, lattice::Lattice, orblocat::AbstractVector{<:ReducedCoordinates})
    size(A) == (numorb(hrh), numorb(hrh), 3) || error("Buffer size mismatch.")
    length(orblocat) == numorb(hrh) || error("Mismatched orbital locations.")
    _Arh_calc!(A, hrh.core, hrh.ij, Val(:partial), k, lattice, orblocat)
    for i in 1:3
        _hrh_fill!(view(A, :, :, i), numorb(hrh), Val(hrh.uplo))
    end
    return A
end

@inline function _hrh_fill!(A, norb, ::Val{:U})
    @inbounds for j in 1:(norb-1), i in (j+1):norb
        A[i, j] = conj(A[j, i])
    end
    @inbounds for i in 1:norb
        A[i, i] = complex(real(A[i, i]), 0)
    end
end
@inline function _hrh_fill!(A, norb, ::Val{:L})
    @inbounds for j in 2:norb, i in 1:(j-1)
        A[i, j] = conj(A[j, i])
    end
    @inbounds for i in 1:norb
        A[i, i] = complex(real(A[i, i]), 0)
    end
end
function Base.show(io::IO, hrh::HermitianReciprocalHoppings)
    print(io, "HermitianReciprocalHoppings with $(sum(Nhop(hrh))) hoppings and $(numorb(hrh)) orbitals.")
end
Base.eltype(::Type{HermitianReciprocalHoppings{T,U}}) where {T,U} = T
function Base.filter(F, hrh::HermitianReciprocalHoppings)
    brh = filter(F, hrh.core)
    return HermitianReciprocalHoppings(brh, hrh.uplo)
end
function prune(hrh::HermitianReciprocalHoppings, cutoff::Real)
    brh = prune(hrh.core, cutoff)
    return HermitianReciprocalHoppings(brh, hrh.uplo)
end
function Base.real(hrh::HermitianReciprocalHoppings{T,U}) where {T<:Complex,U}
    brh = real(hrh.core)
    return HermitianReciprocalHoppings(brh, hrh.uplo)
end
Base.zero(::Type{HermitianReciprocalHoppings{T,U}}, norb::Integer) where {T,U} =
    HermitianReciprocalHoppings(zero(U, norb), :U)
Base.convert(::Type{HermitianReciprocalHoppings{T₁,U₁}}, hrh::HermitianReciprocalHoppings{T₂,U₂}) where {T₁,U₁,T₂,U₂} =
    HermitianReciprocalHoppings(convert(U₁, hrh.core), hrh.uplo)
"""
	translate(hrh::HermitianReciprocalHoppings, centres::Pair{<:Integer, <:AbstractVector{<:Integer}}...)

You only need to input the wannier centres R which isn't in the [0, 0, 0] cell.

See also `translate!`.
"""
function translate(hrh::HermitianReciprocalHoppings, aim_centres::Pair{<:Integer,<:AbstractVector{<:Integer}}...)
    brh = translate(hrh.core, aim_centres)
    return HermitianReciprocalHoppings(brh, hrh.uplo)
end
Base.merge(hrh::HermitianReciprocalHoppings, Arh::AbstractReciprocalHoppings) =
    HermitianReciprocalHoppings(merge(hrh.core, BaseReciprocalHoppings(Arh)), hrh.uplo)
Base.merge(hrh::HermitianReciprocalHoppings, hr::AbstractHR) =
    HermitianReciprocalHoppings(merge(hrh.core, hr), hrh.uplo)
Base.merge(hrh::HermitianReciprocalHoppings, arrays::Union{AbstractArray{<:Hopping},Hopping}...) =
    HermitianReciprocalHoppings(merge(hrh.core, arrays...), hrh.uplo)
Base.union(hrh::HermitianReciprocalHoppings, Arhs::AbstractReciprocalHoppings...) =
    HermitianReciprocalHoppings(union(BaseReciprocalHoppings(hrh), BaseReciprocalHoppings.(Arhs)...), hrh.uplo)

function spinorize(hrh::HermitianReciprocalHoppings{T,U}; mode=conj) where {T,U}
    return HermitianReciprocalHoppings(spinorize(hrh.core; mode), hrh.uplo)
end

function HermitianReciprocalHoppings(hr::HRFull{T}, uplo::Symbol=:U) where {T}
    ij = _triangular_index(uplo, numorb(hr))
    brh = BaseReciprocalHoppings(hr)
    return HermitianReciprocalHoppings{T,BaseReciprocalHoppings{T}}(brh, uplo, ij)
end
function HermitianReciprocalHoppings(hr::HRHermitian{T}, uplo::Symbol=:U) where {T}
    ij = _triangular_index(:U, numorb(hr))
    brh = BaseReciprocalHoppings(hr)
    return HermitianReciprocalHoppings{T,BaseReciprocalHoppings{T}}(brh, :U, ij)
end
function HermitianReciprocalHoppings(pair_int::PairIntHRFull{T}, uplo::Symbol=:U) where {T}
    ij = _triangular_index(uplo, pair_int.np)
    brh = BaseReciprocalHoppings(pair_int)
    return HermitianReciprocalHoppings{T,BaseReciprocalHoppings{T}}(brh, uplo, ij)
end
function HermitianReciprocalHoppings(pair_int::PairIntHRHermitian{T}, uplo::Symbol=:U) where {T}
    ij = _triangular_index(:U, pair_int.np)
    brh = BaseReciprocalHoppings(pair_int)
    return HermitianReciprocalHoppings{T,BaseReciprocalHoppings{T}}(brh, :U, ij)
end
function HermitianReciprocalHoppings(brh::BaseReciprocalHoppings{T}, uplo::Symbol=:U) where {T}
    ij = _triangular_index(uplo, numorb(brh))
    return HermitianReciprocalHoppings{T,BaseReciprocalHoppings{T}}(brh, uplo, ij)
end
function HermitianReciprocalHoppings(rh::ReciprocalHoppings{T,U}, uplo::Symbol=:U) where {T,U}
    ij = _triangular_index(uplo, numorb(rh))
    brh = BaseReciprocalHoppings(rh)
    return HermitianReciprocalHoppings{T,U}(brh, uplo, ij)
end
# function HermitianReciprocalHoppings(hrh::HermitianReciprocalHoppings{T,U}, uplo::Symbol=:U) where {T,U}
#     uplo == hrh.uplo && return hrh
#     ij = _triangular_index(uplo, numorb(hrh))
#     return HermitianReciprocalHoppings{T,U}(hermitian(hrh.core), uplo, ij)
# end
function _triangular_index(uplo::Symbol, norb)
    if uplo ∉ (:U, :L)
        throw(ArgumentError("uplo must be either :U or :L, got :$uplo"))
    end
    return _triangular_index(Val(uplo), norb)
end
function _triangular_index(::Val{:U}, norb)
    n = 0
    ij = Vector{Tuple{Int,Int}}(undef, norb * (norb + 1) ÷ 2)
    @inbounds for j in 1:norb, i in 1:j
        n += 1
        ij[n] = (i, j)
    end
    return ij
end
function _triangular_index(::Val{:L}, norb)
    n = 0
    ij = Vector{Tuple{Int,Int}}(undef, norb * (norb + 1) ÷ 2)
    @inbounds for j in 1:norb, i in j:norb
        n += 1
        ij[n] = (i, j)
    end
    return ij
end
function HR(hrh::HermitianReciprocalHoppings{T,U}) where {T,U}
    N = sum(Nhop(hrh))
    is = Vector{Int}(undef, N)
    js = Vector{Int}(undef, N)
    Rs = Vector{ReducedCoordinates{Int}}(undef, N)
    vals = Vector{T}(undef, N)

    norb = numorb(hrh)
    n = 0
    if hrh.uplo == :U
        for j in 2:norb, i in 1:(j-1)
            for hop in hops(hrh)[i, j]
                n += 1
                is[n] = i
                js[n] = j
                Rs[n] = hop.R
                vals[n] = hop.t
            end
        end
    else
        for j in 1:(norb-1), i in (j+1):norb
            for hop in hops(hrh)[i, j]
                n += 1
                is[n] = j
                js[n] = i
                Rs[n] = -hop.R
                vals[n] = conj(hop.t)
            end
        end
    end
    for i in 1:norb
        for hop in hops(hrh)[i, i]
            _isnegative_vector(R[iR]) && continue
            n += 1
            is[n] = i
            js[n] = i
            Rs[n] = hop.R
            vals[n] = hop.t
        end
    end
    # ip==jp的项存在多余的，N的计数偏大。
    resize!(is, n)
    resize!(js, n)
    resize!(Rs, n)
    resize!(vals, n)

    return HRHermitian(is, js, Rs, vals, norb; inbound=true, dofilter=false)
end
