function is_diagonal(A)::Bool
	(m, n) = size(A)
	(m == n) || return false  # 非方阵直接不是对角阵
	@inbounds for j in 1:n, i in 1:n
		i ≠ j && !iszero(A[i, j]) && return false
	end
	return true
end
