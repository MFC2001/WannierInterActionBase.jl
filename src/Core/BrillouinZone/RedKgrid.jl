"""
	RedKgrid <: AbstractBrillouinZone

Actually it's a simple expansion of [`MonkhorstPack`](@ref).
It supports bracket-based access.

```julia
julia> kgrid = MonkhorstPack([5, 5, 1]; kshift = [1//2, 0, 0])
MonkhorstPack([5, 5, 1], [0.5, 0.0, 0.0])

julia> kgrid = RedKgrid(kgrid)
Reducibal kgrids with 25 reducible k-points.

julia> length(kgrid)
25

julia> kgrid[2]
3-element ReducedCoordinates{Rational{Int64}} with indices SOneTo(3):
 -3//10
 -2//5
   0

```
"""
struct RedKgrid <: AbstractBrillouinZone
	kdirect::Vector{ReducedCoordinates{Rational{Int}}}
	kgrid_size::Vec3{Int}
	kshift::Vec3{Rational{Int}}
end
function RedKgrid(kgrid::MonkhorstPack)
	start = -floor.(Int, (kgrid.kgrid_size .- 1) .// 2)
	stop = ceil.(Int, (kgrid.kgrid_size .- 1) .// 2)

	kgrid_index = [-kgrid.kshift .+ Vec3([i, j, k]) for i ∈ start[1]:stop[1], j ∈ start[2]:stop[2], k ∈ start[3]:stop[3]]
	kgrid_index = reshape(kgrid_index, :)
	kdirect = [x .// kgrid.kgrid_size for x in kgrid_index]

	return RedKgrid(normalize_kdirect.(kdirect), kgrid.kgrid_size, kgrid.kshift)
end
function Base.show(io::IO, redkgrid::RedKgrid)
	print(io, "Reducible kgrids with $(length(redkgrid.kdirect)) reducible k-points.")
end
Base.iterate(redkgrid::RedKgrid, state = 1) = state > length(redkgrid) ? nothing : (redkgrid[state], state + 1)
Base.eltype(::RedKgrid) = ReducedCoordinates{Rational{Int}}
Base.length(redkgrid::RedKgrid) = length(redkgrid.kdirect)
# Base.size(redkgrid::RedKgrid) = size(redkgrid.kdirect)
Base.getindex(redkgrid::RedKgrid, index...) = getindex(redkgrid.kdirect, index...)
Base.setindex!(redkgrid::RedKgrid, v, i::Int) = (redkgrid.kdirect[i] = v)
Base.firstindex(::RedKgrid) = 1
Base.lastindex(redkgrid::RedKgrid) = length(redkgrid)
Base.eachindex(redkgrid::RedKgrid) = eachindex(redkgrid.kdirect)
Base.keys(redkgrid::RedKgrid) = keys(redkgrid.kdirect)
