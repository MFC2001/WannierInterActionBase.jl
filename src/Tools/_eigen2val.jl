function _eigen2val(eigens::AbstractArray{T}) where {T <: Eigen}
	isempty(eigens) && return Matrix{Int}(undef, 0, 0)
	len_values = [length(e.values) for e in eigens]
	max_len = maximum(len_values)
	n = length(eigens)
	value = fill(convert(eltype(first(eigens).values), NaN), max_len, n)
	# Matrix{eltype(first(eigens).values)}(undef, max_len, n)
	for (i, e) in enumerate(eigens)
		value[1:len_values[i], i] .= e.values
	end
	return value
end
function _eigen2val(eigens::T) where {T <: Eigen}
	return eigens.values
end