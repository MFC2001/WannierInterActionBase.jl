
# useful for reduced coordinates.
_is_approx_integer(r; atol = 1e-6) = all(ri -> abs(ri - round(ri)) ≤ atol, r)

# function calc_vecangle(v1::AbstractVector{<:Real}, v2::AbstractVector{<:Real})::Real
# 	return atan(norm(v1 × v2), v1 ⋅ v2)
# end
@inline norm2(v) = v ⋅ v

# used for print(::Hopping)
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
