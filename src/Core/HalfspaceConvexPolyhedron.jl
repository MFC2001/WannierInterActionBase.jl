export HalfspaceConvexPolyhedron, is_inside, is_interior, is_inside_detail, is_interior_detail
"""
	HalfspaceConvexPolyhedron(normals::AbstractVector, offsets::AbstractVector) -> HalfspaceConvexPolyhedron

Can be used to define FBZ or miniFBZ.
"""
struct HalfspaceConvexPolyhedron{T <: Real, V <: AbstractVector{T}}
	normals::Vector{V}
	offsets::Vector{T}
	function HalfspaceConvexPolyhedron(normals::Vector{V}, offsets::Vector{T}) where {T <: Real, V <: AbstractVector{T}}
		length(normals) == length(offsets) || throw(DimensionMismatch("The number of normals must match the number of offsets"))
		any(iszero, normals) && throw(ArgumentError("normal shouldn't be zero vector"))
		new{T, V}(normals, offsets)
	end
end
HalfspaceConvexPolyhedron(normals, offsets) = HalfspaceConvexPolyhedron(collect(normals), collect(offsets))

function point_in_halfspaces(
	p::AbstractVector,
	poly::HalfspaceConvexPolyhedron;
	open_set::Bool = false,
	tol::Real = sqrt(eps(Float64)),
)::Bool
	if open_set
		for (n, b) in zip(poly.normals, poly.offsets)
			val = dot(p, n)
			val ≥ b - tol && return false
		end
	else
		for (n, b) in zip(poly.normals, poly.offsets)
			val = dot(p, n)
			val > b + tol && return false
		end
	end
	return true
end
function point_in_halfspaces_detail(
	p::AbstractVector,
	poly::HalfspaceConvexPolyhedron{T, V};
	open_set::Bool = false,
	tol::Real = sqrt(eps(Float64)),
)::Tuple{Bool, V} where {T, V}
	if open_set
		for (n, b) in zip(poly.normals, poly.offsets)
			val = dot(p, n)
			val ≥ b - tol && return (false, n)
		end
	else
		for (n, b) in zip(poly.normals, poly.offsets)
			val = dot(p, n)
			val > b + tol && return (false, n)
		end
	end
	center = similar(poly.normals[begin])
	fill!(center, zero(T))
	return (true, center) # only used for FBZ.
end

Base.in(p::AbstractVector, poly::HalfspaceConvexPolyhedron) = point_in_halfspaces(p, poly)
is_inside(p::AbstractVector, poly::HalfspaceConvexPolyhedron; kwargs...) =
	point_in_halfspaces(p, poly; kwargs...)
is_interior(p::AbstractVector, poly::HalfspaceConvexPolyhedron; tol = eps(Float64)) =
	point_in_halfspaces(p, poly; open_set = true, tol)
is_inside_detail(p::AbstractVector, poly::HalfspaceConvexPolyhedron; kwargs...) =
	point_in_halfspaces_detail(p, poly; kwargs...)
is_interior_detail(p::AbstractVector, poly::HalfspaceConvexPolyhedron; tol = eps(Float64)) =
	point_in_halfspaces_detail(p, poly; open_set = true, tol)
