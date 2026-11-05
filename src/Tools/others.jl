
# useful for reduced coordinates.
"""
	_is_approx_integer(r; atol=1e-6) -> Bool

Return `true` if every component of `r` is an integer within `atol`.
"""
_is_approx_integer(r; atol = 1e-6) = all(ri -> abs(ri - round(ri)) ≤ atol, r)

# function calc_vecangle(v1::AbstractVector{<:Real}, v2::AbstractVector{<:Real})::Real
# 	return atan(norm(v1 × v2), v1 ⋅ v2)
# end
"""
	norm2(v) -> Real

Return the squared Euclidean norm `v ⋅ v` of `v`.
"""
@inline norm2(v) = v ⋅ v

"""
	subscriptnumber(i::Int) -> String

Return `i` written with Unicode subscript digits, with a leading minus sign for negative
numbers.
"""
function subscriptnumber(i::Int)
	if i < 0
		c = [Char(0x208B)]
	else
		c = []
	end
	for d in reverse(digits(abs(i)))
		push!(c, Char(0x2080 + d))
	end
	return join(c)
end

"""
	superscriptnumber(i::Int) -> String

Return `i` written with Unicode superscript digits, with a leading minus sign for
negative numbers.
"""
function superscriptnumber(i::Int)
	if i < 0
		c = [Char(0x207B)]
	else
		c = []
	end
	for d in reverse(digits(abs(i)))
		if d == 0
			push!(c, Char(0x2070))
		end
		if d == 1
			push!(c, Char(0x00B9))
		end
		if d == 2
			push!(c, Char(0x00B2))
		end
		if d == 3
			push!(c, Char(0x00B3))
		end
		if d > 3
			push!(c, Char(0x2070 + d))
		end
	end
	return join(c)
end
