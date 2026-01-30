
abstract type AbstractRealSpaceWaveFunction end


struct wannier{T <: Number} <: AbstractRealSpaceWaveFunction
    center
	origin::Vec3{Float64}
	grid_vecs::Lattice{Float64}
	grid_size::Vec3{Int}
	value::Array{T, 3}
end
