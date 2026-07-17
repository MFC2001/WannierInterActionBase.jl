function _kline_shiftΓ!(kline::Kline; atol::Real = 1e-6, rtol::Real = 1e-5)
	ΓI = findall(iszero, kline)
	isempty(ΓI) && kline

	n = length(kline)
	for i in ΓI
		if i == 1
			next = kline[2]
			kline[1] = next .* min(atol / norm(next), rtol)
			continue
		elseif i == n
			prev = kline[n-1]
			kline[n] = prev .* min(atol / norm(prev), rtol)
			continue
		end

		prev = kline[i-1]
		next = kline[i+1]

		# 判断三点是否共线（两侧点位于同一条线段上）
		cross_prod = cross(prev, next)
		is_collinear = all(c -> abs(c) < 1e-10, cross_prod)

		if is_collinear
			kline[i] = next .* min(atol / norm(next), rtol)
		else
			middle = (prev + next) ./ 2
			kline[i] = middle .* min(atol / norm(middle), rtol)
		end
	end

	return kline
end
