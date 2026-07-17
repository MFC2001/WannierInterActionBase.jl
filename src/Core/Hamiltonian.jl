
abstract type AbstractHamilton end

struct TBHamiltonian{T <: Real, HT <: Number} <: AbstractHamilton
	kpoint::Vec3{T}
	H::Matrix{HT}
	ittr::Vector{CartesianIndex{2}}
	uplo::Symbol
end

Base.getindex(H::TBHamiltonian,kcoord::AbstractVector{<:Real}) = Hamiltonian(kcoord,)

function Hamiltonian(kpoint::AbstractVector{<:Real}, Hk::AbstractMatrix{<:Number}, uplo = :U; ittr = nothing)::Hamiltonian
	if isnothing(ittr)
		ittr = Vector{CartesianIndex{2}}(undef, Int(norb * (norb + 1) / 2))
		n = 0
		if uplo == :U
			for j in axes(Hk, 1), i in Base.OneTo(j)
				n += 1
				ittr[n] = CartesianIndex(i, j)
			end
		elseif uplo == :L
			for i in axes(Hk, 1), j in i:norb
				n += 1
				ittr[n] = CartesianIndex(i, j)
			end
		end
	end

	return Hamiltonian(kpoint, Hk, ittr, uplo)
end

Base.:(+)(H₁::Hamiltonian, H₂::Hamiltonian) = 1

eigen(Hk::Hamiltonian; permute::Bool = true, scale::Bool = true) = eigen(Hermitian(Hk.H, Hk.uplo); permute, scale)
