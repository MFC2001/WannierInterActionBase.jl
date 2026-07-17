
struct DependentArray{T, N} <: AbstractArray{T, N} where {T, N}
	mapping::Function
	independent::Vector{T}
	function DependentArray{T}(mapping::Function, independent::AbstractVector{T}) where {T}
		N = length(first(methods(mapping)).sig.parameter) - 2
		return new{T, N}(mapping::Function, independent::AbstractVector{T})
	end
end

const DependentVector{T} = DependentArray{T, 1}
const DependentMatrix{T} = DependentArray{T, 2}




Base.length(A::DependentArray{T, N}) where {T, N} = length(A.independent)
Base.size(A::DependentArray{T, N}) where {T, N} = d <= N ? size(A)[d] : 1
Base.size(A::DependentArray{T, N}, d::Integer) where {T, N} = d <= N ? size(A)[d] : 1
Base.getindex(A::DependentArray{T, N}, i::Integer)::T where {T, N} = A.independent[i]
Base.getindex(A::DependentArray{T, N}, i::Integer, Is::Integer...)::T where {T, N} = A.mapping(A, i, Is...)

