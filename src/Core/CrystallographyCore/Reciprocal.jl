
export ReciprocalLattice, reciprocal, isreciprocal

"""
	ReciprocalLattice(data::AbstractMatrix)

Construct a `ReciprocalLattice` from a matrix.

!!! note
	The basis vectors of the matrix are stored as columns.

!!! warning
	Avoid using this constructor directly. Use `reciprocal` instead.
"""
@struct_hash_equal_isequal_isapprox struct ReciprocalLattice{T} <: AbstractLattice{T}
	data::Mat3{T}
end
ReciprocalLattice(data::AbstractMatrix) = ReciprocalLattice(Mat3(data))
ReciprocalLattice(𝐚::AbstractVector, 𝐛::AbstractVector, 𝐜::AbstractVector) = ReciprocalLattice(hcat(𝐚, 𝐛, 𝐜))

"""
	basisvectors(lattice::ReciprocalLattice)

Get the three basis vectors from a reciprocal lattice.
"""
basisvectors(lattice::ReciprocalLattice) =
	lattice[begin:(begin+2)], lattice[(begin+3):(begin+5)], lattice[(begin+6):end]

"""
	reciprocal(lattice::Lattice)
	reciprocal(lattice::ReciprocalLattice)

Get the reciprocal of a `Lattice` or a `ReciprocalLattice`.
"""
function reciprocal(lattice::Lattice)
	Ω = det(lattice.data)  # Cannot use `cellvolume`, it takes the absolute value!
	𝐚, 𝐛, 𝐜 = basisvectors(lattice)
	return 2π * inv(Ω) * ReciprocalLattice(hcat(cross(𝐛, 𝐜), cross(𝐜, 𝐚), cross(𝐚, 𝐛)))
end
function reciprocal(lattice::ReciprocalLattice)
	Ω⁻¹ = det(lattice.data)  # Cannot use `cellvolume`, it takes the absolute value!
	𝐚⁻¹, 𝐛⁻¹, 𝐜⁻¹ = basisvectors(lattice)
	return 2π * inv(Ω⁻¹) * Lattice(hcat(cross(𝐛⁻¹, 𝐜⁻¹), cross(𝐜⁻¹, 𝐚⁻¹), cross(𝐚⁻¹, 𝐛⁻¹)))
end

isreciprocal(a::ReciprocalLattice, b::Lattice) = parent(a)' * parent(b) ≈ I
isreciprocal(a::Lattice, b::ReciprocalLattice) = isreciprocal(b, a)

# See https://github.com/JuliaLang/julia/blob/v1.10.0-beta1/stdlib/LinearAlgebra/src/uniformscaling.jl#L130-L131
Base.one(::Type{ReciprocalLattice{T}}) where {T} =
	ReciprocalLattice(Mat3(SDiagonal(one(T), one(T), one(T))))
Base.one(lattice::ReciprocalLattice) = one(typeof(lattice))

# See https://github.com/JuliaLang/julia/blob/v1.10.0-beta1/stdlib/LinearAlgebra/src/uniformscaling.jl#L132-L133
Base.oneunit(::Type{ReciprocalLattice{T}}) where {T} =
	ReciprocalLattice(Mat3(SDiagonal(oneunit(T), oneunit(T), oneunit(T))))
Base.oneunit(lattice::ReciprocalLattice) = oneunit(typeof(lattice))

# Similar to https://github.com/JuliaCollections/IterTools.jl/blob/0ecaa88/src/IterTools.jl#L1028-L1032
Base.iterate(iter::ReciprocalLattice, state = 1) = iterate(parent(iter), state)

Base.IteratorSize(::Type{<:ReciprocalLattice}) = Base.HasShape{2}()

Base.eltype(::Type{ReciprocalLattice{T}}) where {T} = T

Base.length(::ReciprocalLattice) = 9

Base.size(::ReciprocalLattice) = (3, 3)
# See https://github.com/rafaqz/DimensionalData.jl/blob/bd28d08/src/array/array.jl#L74
Base.size(lattice::ReciprocalLattice, dim) = size(parent(lattice), dim)  # Here, `parent(A)` is necessary to avoid `StackOverflowError`.

Base.parent(lattice::ReciprocalLattice) = lattice.data

Base.getindex(lattice::ReciprocalLattice, i...) = getindex(parent(lattice), i...)

Base.firstindex(::ReciprocalLattice) = 1

Base.lastindex(::ReciprocalLattice) = 9

# You need this to let the broadcasting work.
Base.:*(lattice::ReciprocalLattice, reduced::ReducedCoordinates) = CartesianCoordinates(parent(lattice) * reduced)

# You need this to let the broadcasting work.
Base.:\(lattice::ReciprocalLattice, cartesian::CartesianCoordinates) = ReducedCoordinates(parent(lattice) \ cartesian)

# You need this to let the broadcasting work.
Base.:*(lattice::ReciprocalLattice, x::Number) = ReciprocalLattice(parent(lattice) * x)
Base.:*(x::Number, lattice::ReciprocalLattice) = lattice * x

# You need this to let the broadcasting work.
Base.:/(lattice::ReciprocalLattice, x::Number) = ReciprocalLattice(parent(lattice) / x)
Base.:/(::Number, ::ReciprocalLattice) =
	throw(ArgumentError("you cannot divide a number by a reciprocal lattice!"))

Base.:+(lattice::ReciprocalLattice) = lattice

Base.:-(lattice::ReciprocalLattice) = ReciprocalLattice(-parent(lattice))

Base.convert(::Type{ReciprocalLattice{T}}, lattice::ReciprocalLattice{T}) where {T} =
	lattice
Base.convert(::Type{ReciprocalLattice{T}}, lattice::ReciprocalLattice{S}) where {S, T} =
	ReciprocalLattice(convert(Mat3{T}, parent(lattice)))

# See https://github.com/JuliaLang/julia/blob/v1.10.0-beta3/base/refpointer.jl#L95-L96
Base.ndims(::Type{<:ReciprocalLattice}) = 2
Base.ndims(::ReciprocalLattice) = 2

# See https://github.com/JuliaLang/julia/blob/v1.10.0-rc2/base/broadcast.jl#L741
Base.broadcastable(lattice::ReciprocalLattice) = lattice

# See https://github.com/JuliaLang/julia/blob/v1.10.0-rc2/base/broadcast.jl#L49
Base.BroadcastStyle(::Type{<:ReciprocalLattice}) = Broadcast.Style{ReciprocalLattice}()

# See https://github.com/JuliaLang/julia/blob/v1.10.0-rc2/base/broadcast.jl#L135
Base.BroadcastStyle(
	::Broadcast.AbstractArrayStyle{0}, b::Broadcast.Style{ReciprocalLattice},
) = b

# See https://github.com/JuliaLang/julia/blob/v1.10.0-rc2/base/broadcast.jl#L1114-L1119
Base.copy(bc::Broadcast.Broadcasted{Broadcast.Style{ReciprocalLattice}}) =
	ReciprocalLattice(Mat3(x for x in bc))  # For uniary and binary functions

Base.broadcasted(::typeof(/), ::Number, ::ReciprocalLattice) =
	throw(ArgumentError("you cannot divide a number by a reciprocal lattice!"))
