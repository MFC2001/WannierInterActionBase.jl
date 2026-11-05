"""
	plane_frame(normal; reference=[0, 0, 1], origin=[0, 0, 0], atol=1e-8)

Return the orthonormal frame `(O=origin, ê₁, ê₂, ê₃)` of the plane with the given
`normal`: `ê₃` is parallel to `normal` and `ê₁` is built from `reference`. If
`reference` is parallel to `normal`, a coordinate axis is used instead.
"""
function plane_frame(normal::AbstractVector{<:Real};
    reference::AbstractVector{<:Real}=[0, 0, 1], origin::AbstractVector{<:Real}=[0, 0, 0],
    atol::Real=1e-8,
)
    normal = Vec3(normal)
    normal = normal ./ norm(normal)
    reference = Vec3(reference)
    reference = reference ./ norm(reference)
    origin = Vec3(origin)

    ê₃ = normal
    if isapprox(ê₃ ⋅ reference, 1; atol)
        (_, i) = findmin(abs, ê₃)
        reference = [0, 0, 0]
        reference[i] = 1
        reference = Vec3(reference)
    end

    ê₁ = reference - (reference ⋅ ê₃) * ê₃
    ê₁ = ê₁ ./ norm(ê₁)
    ê₂ = ê₃ × ê₁

    return (O=origin, ê₁, ê₂, ê₃)
end