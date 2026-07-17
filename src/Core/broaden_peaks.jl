export AbstractBroadeningKernel, broaden_peaks
abstract type AbstractBroadeningKernel end
function broaden_peaks(sym::Symbol, x::AbstractVector{<:Real},
	x₀::AbstractVector{<:Real},
	A::AbstractVector{<:Real} = ones(Float64, length(x₀));
	normalize::Bool = true,
	scale::Real = 0.1, scale_vector::AbstractVector{<:Real} = fill(scale, length(x₀)),
)
	x = collect(x)
	y = zeros(Float64, length(x))
	@assert length(x₀) == length(A) == length(scale_vector) "Mismatched length!"
	for (xx₀, AA, ss) in zip(x₀, A, scale_vector)
		f = BroadeningKernel(sym, xx₀, ss, AA; normalize)
		y .+= f.(x)
	end
	return y
end
export BroadeningKernel, Lorentzian, Gaussian
BroadeningKernel(sym::Symbol, args...; kwargs...) = BroadeningKernel(Val(sym), args...; kwargs...)
"""
	Lorentzian(x₀::Real, γ::Real; normalize::Bool = true)
	Lorentzian(x₀::Real, γ::Real, A::Real; normalize::Bool = true)

```julia
h = normalize ? A / (π * γ) : A
```

```math
L(x) = h \\frac{\\gamma^2}{(x - x_0)^2 + \\gamma^2}
"""
struct Lorentzian <: AbstractBroadeningKernel
	x₀::Float64
	γ::Float64
	h::Float64
	γsquare::Float64
	hγsquare::Float64
end
BroadeningKernel(::Val{:lorentzian}, args...; kwargs...) = Lorentzian(args...; kwargs...)
Lorentzian(x₀::Real, γ::Real; kwargs...) = Lorentzian(x₀, γ, 1.0; kwargs...)
function Lorentzian(x₀::Real, γ::Real, A::Real; normalize::Bool = true)
	@assert γ > 0 "γ must be positive"
	γsquare = γ * γ
	h = normalize ? A / (π * γ) : A
	# if normalize
	# 	∫L(x) dx = A
	# else
	# 	L(x₀) = A
	# end
	return Lorentzian(x₀, γ, h, γsquare, h * γsquare)
end
function (L::Lorentzian)(x::Real)
	dx = x - L.x₀
	return L.hγsquare / (dx * dx + L.γsquare)
end
"""
	Gaussian(x₀::Real, σ::Real; normalize::Bool = true)
	Gaussian(x₀::Real, σ::Real, A::Real; normalize::Bool = true)

```julia
h = normalize ? A / (σ * sqrt(2π)) : A
```

```math
G(x) = h exp{-\\frac{(x - x_0)^2}{2\\sigma^2}}
"""
struct Gaussian <: AbstractBroadeningKernel
	x₀::Float64
	σ::Float64
	h::Float64
	σsquare::Float64
	inv_2σsquare::Float64
end
BroadeningKernel(::Val{:gauss}, args...; kwargs...) = Gaussian(args...; kwargs...)
Gaussian(x₀::Real, σ::Real; kwargs...) = Gaussian(x₀, σ, 1.0; kwargs...)
function Gaussian(x₀::Real, σ::Real, A::Real; normalize::Bool = true)
	@assert σ > 0 "σ must be positive"
	σsquare = σ * σ
	inv_2σsquare = 1 / (2σsquare)   # 预乘系数，调用时避免除法
	h = normalize ? A / (σ * sqrt(2π)) : A
	# if normalize
	#     ∫G(x) dx = A
	# else
	#     G(x₀) = A
	# end
	return Gaussian(x₀, σ, h, σsquare, inv_2σsquare)
end
function (g::Gaussian)(x::Real)
	dx = x - g.x₀
	return g.h * exp(-dx * dx * g.inv_2σsquare)
end

