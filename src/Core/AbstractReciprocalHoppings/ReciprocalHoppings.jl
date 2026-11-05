"""
	ReciprocalHoppings{T, U <: BaseReciprocalHoppings{T}} <: AbstractReciprocalHoppings{T}

```julia
julia> hr isa AbstractHR
true

julia> rh = ReciprocalHoppings(hr);

julia> k isa ReducedCoordinates
true

julia> rh(k);
```
"""
struct ReciprocalHoppings{T,U<:BaseReciprocalHoppings{T}} <: AbstractReciprocalHoppings{T}
    core::U
    ij::Vector{Tuple{Int,Int}}
end

@inline function (rh::ReciprocalHoppings)(k::ReducedCoordinates)
    A = Matrix{ComplexF64}(undef, numorb(rh), numorb(rh))
    _Arh_calc!(A, rh.core, rh.ij, k)
    return A
end
@inline function (rh::ReciprocalHoppings)(A::AbstractMatrix, k::ReducedCoordinates)
    size(A) == (numorb(rh), numorb(rh)) || error("Buffer size mismatch.")
    _Arh_calc!(A, rh.core, rh.ij, k)
    return A
end

@inline function (rh::ReciprocalHoppings)(::Val{:partial}, k::ReducedCoordinates, lattice::Lattice)
    A = Array{ComplexF64}(undef, numorb(rh), numorb(rh), 3)
    _Arh_calc!(A, rh.core, rh.ij, Val(:partial), k, lattice)
    return A
end
@inline function (rh::ReciprocalHoppings)(::Val{:partial}, A::AbstractArray, k::ReducedCoordinates, lattice::Lattice)
    size(A) == (numorb(rh), numorb(rh), 3) || error("Buffer size mismatch.")
    _Arh_calc!(A, rh.core, rh.ij, Val(:partial), k, lattice)
    return A
end

@inline function (rh::ReciprocalHoppings)(k::ReducedCoordinates, orblocat::AbstractVector{<:ReducedCoordinates})
    A = Matrix{ComplexF64}(undef, numorb(rh), numorb(rh))
    length(orblocat) == numorb(rh) || error("Mismatched orbital locations.")
    _Arh_calc!(A, rh.core, rh.ij, k, orblocat)
    return A
end
@inline function (rh::ReciprocalHoppings)(A::AbstractMatrix, k::ReducedCoordinates, orblocat::AbstractVector{<:ReducedCoordinates})
    size(A) == (numorb(rh), numorb(rh)) || error("Buffer size mismatch.")
    length(orblocat) == numorb(rh) || error("Mismatched orbital locations.")
    _Arh_calc!(A, rh.core, rh.ij, k, orblocat)
    return A
end

@inline function (rh::ReciprocalHoppings)(::Val{:partial}, k::ReducedCoordinates, lattice::Lattice, orblocat::AbstractVector{<:ReducedCoordinates})
    A = Array{ComplexF64}(undef, numorb(rh), numorb(rh), 3)
    length(orblocat) == numorb(rh) || error("Mismatched orbital locations.")
    _Arh_calc!(A, rh.core, rh.ij, Val(:partial), k, lattice, orblocat)
    return A
end
@inline function (rh::ReciprocalHoppings)(::Val{:partial}, A::AbstractArray, k::ReducedCoordinates, lattice::Lattice, orblocat::AbstractVector{<:ReducedCoordinates})
    size(A) == (numorb(rh), numorb(rh), 3) || error("Buffer size mismatch.")
    length(orblocat) == numorb(rh) || error("Mismatched orbital locations.")
    _Arh_calc!(A, rh.core, rh.ij, Val(:partial), k, lattice, orblocat)
    return A
end

function Base.show(io::IO, rh::ReciprocalHoppings)
    print(io, "ReciprocalHoppings with $(sum(Nhop(rh))) hoppings and $(numorb(rh)) orbitals.")
end
Base.eltype(::Type{ReciprocalHoppings{T,U}}) where {T,U} = T
function Base.filter(F, rh::ReciprocalHoppings)
    brh = filter(F, rh.core)
    return ReciprocalHoppings(brh)
end
function prune(rh::ReciprocalHoppings, cutoff::Real)
    brh = prune(rh.core, cutoff)
    return ReciprocalHoppings(brh)
end
function Base.real(rh::ReciprocalHoppings{T,U}) where {T<:Complex,U}
    brh = real(rh.core)
    return ReciprocalHoppings(brh)
end
Base.zero(::Type{ReciprocalHoppings{T,U}}, norb::Integer) where {T,U} =
    ReciprocalHoppings(zero(U, norb))
Base.convert(::Type{ReciprocalHoppings{T₁,U₁}}, rh::ReciprocalHoppings{T₂,U₂}) where {T₁,U₁,T₂,U₂} =
    ReciprocalHoppings(convert(U₁, rh.core))
"""
	translate(rh::ReciprocalHoppings, centres::Pair{<:Integer, <:AbstractVector{<:Integer}}...)

You only need to input the wannier centres R which isn't in the [0, 0, 0] cell.

See also `translate!`.
"""
function translate(rh::ReciprocalHoppings, aim_centres::Pair{<:Integer,<:AbstractVector{<:Integer}}...)
    brh = translate(rh.core, aim_centres)
    return ReciprocalHoppings(brh)
end
Base.merge(rh::ReciprocalHoppings, Arh::AbstractReciprocalHoppings) =
    ReciprocalHoppings(merge(rh.core, BaseReciprocalHoppings(Arh)))
Base.merge(rh::ReciprocalHoppings, hr::AbstractHR) =
    ReciprocalHoppings(merge(rh.core, hr))
Base.merge(rh::ReciprocalHoppings, arrays::Union{AbstractArray{<:Hopping},Hopping}...) =
    ReciprocalHoppings(merge(rh.core, arrays...))
Base.union(rh::ReciprocalHoppings, Arhs::AbstractReciprocalHoppings...) =
    ReciprocalHoppings(union(BaseReciprocalHoppings(rh), BaseReciprocalHoppings.(Arhs)...))

function spinorize(rh::ReciprocalHoppings{T,U}; mode=conj) where {T,U}
    return ReciprocalHoppings(spinorize(rh.core; mode))
end

"""
	ReciprocalHoppings(hr::AbstractHR{T}) -> ReciprocalHoppings{T, BaseReciprocalHoppings{T}}
"""
function ReciprocalHoppings(hr::HRFull{T}) where {T}
    brh = BaseReciprocalHoppings(hr)
    ij = _rh_matrix_index(numorb(brh))
    return ReciprocalHoppings{T,BaseReciprocalHoppings{T}}(brh, ij)
end
function ReciprocalHoppings(int::PairIntHRFull{T}) where {T}
    brh = BaseReciprocalHoppings(int)
    ij = _rh_matrix_index(numorb(brh))
    return ReciprocalHoppings{T,BaseReciprocalHoppings{T}}(brh, ij)
end
function ReciprocalHoppings(brh::BaseReciprocalHoppings{T}) where {T}
    ij = _rh_matrix_index(numorb(brh))
    return ReciprocalHoppings{T,BaseReciprocalHoppings{T}}(brh, ij)
end
# function ReciprocalHoppings(Arh::AbstractReciprocalHoppings{T}) where {T}
#     brh = BaseReciprocalHoppings(Arh)
#     ij = _rh_matrix_index(numorb(brh))
#     return ReciprocalHoppings{T,BaseReciprocalHoppings{T}}(brh, ij)
# end
_rh_matrix_index(norb) = reshape([(i, j) for i in 1:norb, j in 1:norb], :)
function HR(rh::ReciprocalHoppings{T,BaseReciprocalHoppings{T}}) where {T}
    N = sum(Nhop(rh))
    is = Vector{Int}(undef, N)
    js = Vector{Int}(undef, N)
    Rs = Vector{ReducedCoordinates{Int}}(undef, N)
    vals = Vector{T}(undef, N)

    norb = numorb(rh)
    n = 0
    for j in 1:norb, i in 1:norb
        for hop in hops(rh)[i, j]
            n += 1
            is[n] = i
            js[n] = j
            Rs[n] = hop.R
            vals[n] = hop.t
        end
    end

    return HRFull(is, js, Rs, vals, norb; inbound=true, dofilter=false)
end
