
"""
	phase_normalize!(mode::Symbol, obj; kwargs...)

Normalize the phases of the eigenvectors in `obj`, in place, and return it. `obj` may be
an `Eigen` object, a vector of `Eigen` objects, a matrix whose columns are vectors, or a
single vector. `mode` selects `:positive_sum` (phase such that the sum is positive real)
or `:positive_maxabs` (phase such that the largest component is positive real).
"""
phase_normalize!(mode::Symbol, band::AbstractVector{<:Eigen}; kwargs...) = map(x -> phase_normalize!(mode, x; kwargs...), band)
"""
	phase_normalize(mode::Symbol, obj; kwargs...)

Return a phase-normalized copy of `obj`; see [`phase_normalize!`](@ref) for the meaning
of `mode` and the accepted types.
"""
phase_normalize(mode::Symbol, band::AbstractVector{<:Eigen}; kwargs...) = map(x -> phase_normalize(mode, x; kwargs...), band)

phase_normalize!(mode::Symbol, band::Eigen; kwargs...) = phase_normalize!(mode, band.vectors; kwargs...)
phase_normalize(mode::Symbol, band::Eigen; kwargs...) = Eigen(band.values, phase_normalize(mode, band.vectors; kwargs...))

function phase_normalize!(mode::Symbol, vectors::AbstractMatrix{<:Number}; kwargs...)
    for i in axes(vectors, 2)
        phase_normalize!(Val(mode), view(vectors, :, i); kwargs...)
    end
    return vectors
end
function phase_normalize(mode::Symbol, vectors::AbstractMatrix{<:Number}; kwargs...)
    vectors = copy(vectors)
    for i in axes(vectors, 2)
        phase_normalize!(Val(mode), view(vectors, :, i); kwargs...)
    end
    return vectors
end
phase_normalize!(mode::Symbol, vector::AbstractVector{<:Number}; kwargs...) = phase_normalize!(Val(mode), vector; kwargs...)
phase_normalize(mode::Symbol, vector::AbstractVector{<:Number}; kwargs...) = phase_normalize!(Val(mode), copy(vector); kwargs...)
function phase_normalize!(::Val{:positive_sum}, vector::AbstractVector{<:Complex}; atol::Real=1e-12)
    isempty(vector) && return vector
    t = sum(vector)
    if abs(imag(t)) > atol # make sure abs(t) is't zero.
        e_neg_iθ = conj(t) / abs(t)
        vector .*= e_neg_iθ
    end
    return vector
end
function phase_normalize!(::Val{:positive_sum}, vector::AbstractVector{<:Real}; atol::Real=1e-12)
    isempty(vector) && return vector
    t = sum(vector)
    if abs(t) > atol && t < 0
        vector .*= -1
    end
    return vector
end
function phase_normalize!(::Val{:positive_maxabs}, vector::AbstractVector{<:Complex}; atol::Real=1e-12)
    isempty(vector) && return vector
    # define maxabs
    k = firstindex(vector)
    maxabs = abs(vector[k])
    for i in eachindex(vector)
        a = abs(vector[i])
        if a > maxabs
            maxabs = a
            k = i
        end
    end
    # process
    if maxabs > atol
        vk = vector[k]
        e_neg_iθ = conj(vk) / maxabs
        vector .*= e_neg_iθ
    end
    return vector
end
function phase_normalize!(::Val{:positive_maxabs}, vector::AbstractVector{<:Real}; atol::Real=1e-12)
    isempty(vector) && return vector
    # define maxabs
    k = firstindex(vector)
    maxabs = abs(vector[k])
    for i in eachindex(vector)
        a = abs(vector[i])
        if a > maxabs
            maxabs = a
            k = i
        end
    end
    if maxabs > atol && vector[k] < 0
        vector .*= -1
    end
    return vector
end
