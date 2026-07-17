"""
	_ispositive_vector(vector; atol::Real = 1e-10)::Bool

Check if the first element of the vector is positive within a tolerance.
"""
function _ispositive_vector(vector; atol::Real = 1e-10)::Bool
	for v in vector
		abs(v) > atol && return v > 0
	end
	return false
end
function _ispositive_vector(vector::AbstractArray{<:Integer}; atol::Real = 1e-10)::Bool
	for v in vector
		v ≠ 0 && return v > 0
	end
	return false
end
"""
	_isnegative_vector(vector; atol::Real = 1e-10)::Bool

Check if the first element of the vector is negative within a tolerance.
"""
function _isnegative_vector(vector; atol::Real = 1e-10)::Bool
	for v in vector
		abs(v) > atol && return v < 0
	end
	return false
end
function _isnegative_vector(vector::AbstractArray{<:Integer}; atol::Real = 1e-10)::Bool
	for v in vector
		v ≠ 0 && return v < 0
	end
	return false
end
