
export Lattice, basisvectors

"""
	AbstractLattice{T}

Represent the real lattices and the reciprocal lattices.
"""
abstract type AbstractLattice{T} end
@struct_hash_equal_isequal_isapprox struct Lattice{T} <: AbstractLattice{T}
	data::Mat3{T}
end
"""
	Lattice(data::AbstractMatrix) -> Lattice

Construct a `Lattice` from a matrix.

!!! note
	The basis vectors of the matrix are stored as columns.

# Examples
```julia
julia> Lattice([
		   1.2 4.5 7.8
		   2.3 5.6 8.9
		   3.4 6.7 9.1
	   ])
Lattice{Float64}
 1.2  4.5  7.8
 2.3  5.6  8.9
 3.4  6.7  9.1
```
"""
Lattice(data::AbstractMatrix) = Lattice(Mat3(data))
"""
	Lattice(𝐚::AbstractVector, 𝐛::AbstractVector, 𝐜::AbstractVector) -> Lattice

Construct a `Lattice` from three basis vectors.

# Examples
```julia
julia> 𝐚, 𝐛, 𝐜 = [1.2, 2.3, 3.4], [4.5, 5.6, 6.7], [7.8, 8.9, 9.10];

julia> Lattice(𝐚, 𝐛, 𝐜)
Lattice{Float64}
 1.2  4.5  7.8
 2.3  5.6  8.9
 3.4  6.7  9.1
```
"""
Lattice(𝐚::AbstractVector, 𝐛::AbstractVector, 𝐜::AbstractVector) = Lattice(hcat(𝐚, 𝐛, 𝐜))
"""
	Lattice(data)

Construct a `Lattice` from, e.g., a vector of three basis vectors.

# Examples
```julia
julia> Lattice([[1.2, 2.3, 3.4], [4.5, 5.6, 6.7], [7.8, 8.9, 9.10]])
Lattice{Float64}
 1.2  4.5  7.8
 2.3  5.6  8.9
 3.4  6.7  9.1

julia> Lattice(((1.1, 2.2, 3.1), (4.4, 5.5, 6.5), (7.3, 8.8, 9.9)))
Lattice{Float64}
 1.1  4.4  7.3
 2.2  5.5  8.8
 3.1  6.5  9.9

julia> Lattice((1.1, 2.2, 3.1, 4.4, 5.5, 6.5, 7.3, 8.8, 9.9))
Lattice{Float64}
 1.1  4.4  7.3
 2.2  5.5  8.8
 3.1  6.5  9.9

julia> Lattice(i * 1.1 for i in 1:9)
Lattice{Float64}
 1.1  4.4  7.700000000000001
 2.2  5.5  8.8
 3.3000000000000003  6.6000000000000005  9.9
```
"""
Lattice(data::NTuple{9}) = Lattice(Mat3(data))
Lattice(data::NTuple{3, NTuple{3}}) = Lattice(mapreduce(collect, hcat, data))
Lattice(data::AbstractVector{<:AbstractVector}) = Lattice(hcat(data[1], data[2], data[3]))  # This is faster
Lattice(iter::Base.Generator) = Lattice(Mat3(iter))
function Lattice(data::Union{AbstractVector, Tuple})
	if length(data) == 9
		return Lattice(Mat3(data))
	elseif length(data) == 3
		return Lattice(mapreduce(collect, hcat, data))
	else
		throw(DimensionMismatch("The length of the tuple must be 3 or 9."))
	end
end

"""
	basisvectors(lattice::Lattice)

Get the three basis vectors from a lattice.
"""
basisvectors(lattice::Lattice) =
	lattice[begin:(begin+2)], lattice[(begin+3):(begin+5)], lattice[(begin+6):end]

# See https://github.com/JuliaLang/julia/blob/v1.10.0-beta1/stdlib/LinearAlgebra/src/uniformscaling.jl#L130-L131
Base.one(::Type{Lattice{T}}) where {T} = Lattice(SDiagonal(one(T), one(T), one(T)))
Base.one(lattice::Lattice) = one(typeof(lattice))

# See https://github.com/JuliaLang/julia/blob/v1.10.0-beta1/stdlib/LinearAlgebra/src/uniformscaling.jl#L132-L133
Base.oneunit(::Type{Lattice{T}}) where {T} =
	Lattice(SDiagonal(oneunit(T), oneunit(T), oneunit(T)))
Base.oneunit(lattice::Lattice) = oneunit(typeof(lattice))

# See https://github.com/JuliaLang/julia/blob/v1.10.0-beta1/stdlib/LinearAlgebra/src/uniformscaling.jl#L134-L135
Base.zero(::Type{Lattice{T}}) where {T} = Lattice(zeros(T, 3, 3))
Base.zero(lattice::Lattice) = zero(typeof(lattice))

# Similar to https://github.com/JuliaCollections/IterTools.jl/blob/0ecaa88/src/IterTools.jl#L1028-L1032
Base.iterate(iter::Lattice, state = 1) = iterate(parent(iter), state)

Base.IteratorSize(::Type{<:Lattice}) = Base.HasShape{2}()

Base.eltype(::Type{Lattice{T}}) where {T} = T

Base.length(::Lattice) = 9

Base.size(::Lattice) = (3, 3)
# See https://github.com/rafaqz/DimensionalData.jl/blob/bd28d08/src/array/array.jl#L74
Base.size(lattice::Lattice, dim) = size(parent(lattice), dim)  # Here, `parent(A)` is necessary to avoid `StackOverflowError`.

Base.parent(lattice::Lattice) = lattice.data

Base.getindex(lattice::Lattice, i...) = getindex(parent(lattice), i...)

Base.firstindex(::Lattice) = 1

Base.lastindex(::Lattice) = 9

# You need this to let the broadcasting work.
Base.:*(lattice::Lattice, reduced::ReducedCoordinates) = CartesianCoordinates(parent(lattice) * reduced)
Base.:*(lattice::Lattice, reduced::AbstractVector) = parent(lattice) * reduced

# You need this to let the broadcasting work.
Base.:\(lattice::Lattice, cartesian::CartesianCoordinates) = ReducedCoordinates(parent(lattice) \ cartesian)
Base.:\(lattice::Lattice, cartesian::AbstractVector) = parent(lattice) \ cartesian

# You need this to let the broadcasting work.
Base.:*(lattice::Lattice, x::Number) = Lattice(parent(lattice) * x)
Base.:*(x::Number, lattice::Lattice) = lattice * x

# You need this to let the broadcasting work.
Base.:/(lattice::Lattice, x::Number) = Lattice(parent(lattice) / x)
Base.:/(::Number, ::Lattice) =
	throw(ArgumentError("you cannot divide a number by a lattice!"))

Base.:+(lattice::Lattice) = lattice
# You need this to let the broadcasting work.
Base.:+(lattice::Lattice, x::Number) = Lattice(parent(lattice) .+ x)
Base.:+(x::Number, lattice::Lattice) = lattice + x

Base.:-(lattice::Lattice) = Lattice(-parent(lattice))
# You need this to let the broadcasting work.
Base.:-(lattice::Lattice, x::Number) = Lattice(parent(lattice) .- x)
Base.:-(x::Number, lattice::Lattice) = -lattice + x

Base.convert(::Type{Lattice{T}}, lattice::Lattice{T}) where {T} = lattice
Base.convert(::Type{Lattice{T}}, lattice::Lattice{S}) where {S, T} =
	Lattice(convert(SMatrix{3, 3, T, 9}, parent(lattice)))

# See https://github.com/JuliaLang/julia/blob/v1.10.0-beta3/base/refpointer.jl#L95-L96
Base.ndims(::Type{<:Lattice}) = 2
Base.ndims(::Lattice) = 2

# See https://github.com/JuliaLang/julia/blob/v1.10.0-rc2/base/broadcast.jl#L741
Base.broadcastable(lattice::Lattice) = lattice

# See https://github.com/JuliaLang/julia/blob/v1.10.0-rc2/base/broadcast.jl#L49
Base.BroadcastStyle(::Type{<:Lattice}) = Broadcast.Style{Lattice}()

# See https://github.com/JuliaLang/julia/blob/v1.10.0-rc2/base/broadcast.jl#L135
Base.BroadcastStyle(::Broadcast.AbstractArrayStyle{0}, b::Broadcast.Style{Lattice}) = b

# See https://github.com/JuliaLang/julia/blob/v1.10.0-rc2/base/broadcast.jl#L1114-L1119
Base.copy(bc::Broadcast.Broadcasted{Broadcast.Style{Lattice}}) = Lattice(x for x in bc)  # For uniary and binary functions

Base.broadcasted(::typeof(/), ::Number, ::Lattice) =
	throw(ArgumentError("you cannot divide a number by a lattice!"))
