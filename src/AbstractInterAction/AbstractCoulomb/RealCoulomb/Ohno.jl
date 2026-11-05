struct RealCoulombOhno{V₀T} <: RealCoulomb
    U₀::Float64
    V₀::V₀T
end
function (v::RealCoulombOhno)(r::AbstractVector{<:Real}; atol::Real=1e-10)::Float64
    r = norm(r)
    return r ≤ atol ? v.U₀ : v.U₀ / sqrt(1 + (v.U₀/v.V₀(r; atol)))
end
function (v::RealCoulombOhno)(r::Real; atol::Real=1e-10)::Float64
    signbit(r) && error("r should be positive")
    return r ≤ atol ? v.U₀ : v.U₀ / sqrt(1 + (v.U₀/v.V₀(r; atol)))
end
"""
	RealCoulomb(:ohno, V₀::V₀T; U₀::Real=1) -> RealCoulombOhno

Calculate Ohno Coulomb potential in real space:

```math
V(r) = \\frac{U₀}{\\sqrt{1 + (\\frac{U₀}{V₀(r)})^2}}
```

unit: r(Å), V(eV), U₀(eV), V₀(eV).

# Example

```julia
julia> V₀ = RealCoulomb(:point; ϵ = 1.5)
julia> V = RealCoulomb(:ohno, V₀; U₀ = 10)
julia> r = [1.0, 0.0, 0.0] # or r = 1.0
julia> V(r)
```
"""
function RealCoulomb(::Val{:ohno}, V₀::V₀T; U₀::Real=1.0) where {V₀T}
    @assert U₀ > 0 "U₀ should be positive"
    return RealCoulombOhno{V₀T}(U₀, V₀)
end
