export ReducedCoordinates, CrystalCoordinates, CartesianCoordinates

"""
Reduced Coordinates.
"""
struct ReducedCoordinates{T <: Real} <: FieldVector{3, T}
	x::T
	y::T
	z::T
end
"""
	CrystalCoordinates = ReducedCoordinates
"""
const CrystalCoordinates = ReducedCoordinates

"""
Cartesian Coordinates.
"""
struct CartesianCoordinates{T <: Real} <: FieldVector{3, T}
	x::T
	y::T
	z::T
end

# See https://juliaarrays.github.io/StaticArrays.jl/dev/pages/api/#StaticArraysCore.FieldVector
StaticArrays.similar_type(::Type{<:ReducedCoordinates}, ::Type{T}, s::Size{(3,)}) where {T <: Real} = ReducedCoordinates{T}
StaticArrays.similar_type(::Type{<:ReducedCoordinates}, ::Type{T}, s::Size{(3,)}) where {T <: Complex} = Vec3{T}
StaticArrays.similar_type(::Type{<:CartesianCoordinates}, ::Type{T}, s::Size{(3,)}) where {T <: Real} = CartesianCoordinates{T}
StaticArrays.similar_type(::Type{<:CartesianCoordinates}, ::Type{T}, s::Size{(3,)}) where {T <: Complex} = Vec3{T}

Base.show(io::IO, c::ReducedCoordinates) = print(io, "R[$(c.x), $(c.y), $(c.z)]")
Base.show(io::IO, c::CartesianCoordinates) = print(io, "C[$(c.x), $(c.y), $(c.z)]")
Base.show(io::IO, ::MIME"text/plain", c::ReducedCoordinates) = show(io, c)
Base.show(io::IO, ::MIME"text/plain", c::CartesianCoordinates) = show(io, c)

# Base.convert(::Type{AbstractVector{T}},X::ReducedCoordinates{S}) where {S, T} = 
# Base.convert(::Type{Lattice{T}}, lattice::Lattice{S}) where {S, T} =
# 	Lattice(convert(SMatrix{3, 3, T, 9}, parent(lattice)))
