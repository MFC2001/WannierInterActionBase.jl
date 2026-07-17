
phase_normalize!(band::AbstractVector{<:Eigen}, mode::Symbol; kwargs...) = map(x -> phase_normalize!(x, mode; kwargs...), band)
phase_normalize(band::AbstractVector{<:Eigen}, mode::Symbol; kwargs...) = map(x -> phase_normalize(x, mode; kwargs...), band)

phase_normalize!(band::Eigen, mode::Symbol; kwargs...) = phase_normalize!(band.vectors, mode; kwargs...)
phase_normalize(band::Eigen, mode::Symbol; kwargs...) = Eigen(band.values, phase_normalize(band.vectors, mode; kwargs...))

function phase_normalize!(vectors::AbstractMatrix{<:Number}, mode::Symbol; kwargs...)
	for i in axes(vectors, 2)
		phase_normalize!(view(vectors, :, i), Val(mode); kwargs...)
	end
	return vectors
end
function phase_normalize(vectors::AbstractMatrix{<:Number}, mode::Symbol; kwargs...)
	vectors = copy(vectors)
	for i in axes(vectors, 2)
		phase_normalize!(view(vectors, :, i), Val(mode); kwargs...)
	end
	return vectors
end
phase_normalize!(vector::AbstractVector{<:Number}, mode::Symbol; kwargs...) = phase_normalize!(vector, Val(mode); kwargs...)
phase_normalize(vector::AbstractVector{<:Number}, mode::Symbol; kwargs...) = phase_normalize!(copy(vector), Val(mode); kwargs...)
function phase_normalize!(vector::AbstractVector{<:Complex}, ::Val{:real_sum}; atol::Real = 1e-12)
	t = sum(vector)
	if abs(imag(t)) > atol #Here this judge make sure abs(t) is't zero.
		e_neg_iθ = conj(t) / abs(t)
		vector .*= e_neg_iθ
	end
	return vector
end
function phase_normalize!(vector::AbstractVector{<:Real}, ::Val{:real_sum}; atol::Real = 1e-12)
	return vector
end
