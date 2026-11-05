struct RealCoulombSoftCore <: RealCoulomb
    ϵ::Float64
    a::Float64
    coff::Float64
    a²::Float64
end
function (v::RealCoulombSoftCore)(r::AbstractVector{<:Real}; atol::Real=1e-10)::Float64
    r² = r ⋅ r
    return v.coff / sqrt(r² + v.a²)
end
function (v::RealCoulombSoftCore)(r::Real; atol::Real=1e-10)::Float64
    signbit(r) && error("r should be positive")
    return v.coff / sqrt(r * r + v.a²)
end
"""
	RealCoulomb(:softcore; ϵ::Real=1, a::Real=0.0) -> RealCoulombSoftCore

Calculate soft-core Coulomb potential in real space:

```math
V(r) = \\frac{e²}{4πϵ₀ϵ} \\frac{1}{\\sqrt{r^2 + a^2}}
```

unit: r(Å), V(eV), ϵ(dimensionless), a(Å).

# Example

```julia
julia> V = RealCoulomb(:softcore; ϵ = 1.5, a = 1)
julia> r = [1.0, 0.0, 0.0] # or r = 1.0
julia> V(r)
```
"""
function RealCoulomb(::Val{:softcore}; ϵ::Real=1.0, a::Real=0.0, atol::Real=1e-10)
    @assert ϵ > 0 "ϵ should be positive"
    @assert a > 0 "ϵ should be positive"
    coff = CoulombScale / ϵ
    if a ≤ atol
        return RealCoulomb(Val(:point); ϵ, a)
    else
        return RealCoulombSoftCore(ϵ, a, coff, a * a)
    end
end
