
#TODO Need more design.


@enum Spin Up = 1 Down = -1 None = 0
Base.:(*)(spin::Spin, i::Integer) = Spin(Int(spin) * i)
Base.:(*)(i::Integer, spin::Spin) = Spin(Int(spin) * i)

mutable struct WannierBasis{T <: Real, VT <: AbstractVector{T}}
	center::VT
	name::String
	σ::Spin
end

nameof(wb::WannierBasis{T}) where {T} = wb.name
center(wb::WannierBasis{T}) where {T} = wb.center
spinof(wb::WannierBasis{T}) where {T} = wb.σ

Base.isequal(wb₁::WannierBasis{T₁}, wb₂::WannierBasis{T₂}) where {T₁, T₂} =
	wb₁.center == wb₂.center && wb₁.name == wb₂.name && wb₁.σ == wb₂.σ
Base.isapprox(wb₁::WannierBasis{T₁}, wb₂::WannierBasis{T₂}; atol = 1e-4) where {T₁, T₂} =
	isapprox(wb₁.center, wb₂.center; atol) && wb₁.name == wb₂.name && wb₁.σ == wb₂.σ
