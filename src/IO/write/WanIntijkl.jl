
function Base.write(io::IO, intvalue::AbstractArray{VT, 4}, ::Type{WanIntijkl0}; comment = "") where {VT <: Number}
	now = Dates.format(Dates.now(), "yyyy/m/d H:M:S")
	comment = comment * " " * now

	println(io, comment)

	nw = size(intvalue, 1)
	@printf(io, "%12u\n", nw)

	for l in 1:nw, k in 1:nw, j in 1:nw, i in 1:nw
		@printf(io, "%7u %7u %7u %7u %12.6f %12.6f\n", i, j, k, l, reim(intvalue[i, j, k, l])...)
	end

	return nothing
end

function Base.write(io::IO, int::WanIntijklR{T}, ::Type{WanIntijklR}; comment = "") where {T <: Number}
	now = Dates.format(Dates.now(), "yyyy/m/d H:M:S")
	comment = comment * " " * now

	println(io, comment)

	@printf(io, "%12u %12u %12u %16u\n", int.nw, int.np, int.nR, int.nv)
	println(io, "pair: ip    iw1    iw2   Rx   Ry   Rz")
	for (ip, pair) in enumerate(int.pair)
		(iw1, iw2, ipR) = pair
		pR = int.pR[ipR]
		@printf(io, "%7u %6u %6u %4u %4u %4u\n", ip, iw1, iw2, pR...)
	end

	println(io, "R:    iR    Rx    Ry    Rz")
	for iR in 1:int.nR
		R = int.R[iR]
		@printf(io, "%7u %5u %5u %5u\n", iR, R...)
	end

	println(io, "int: ip1    ip2    iR    real     imag")
	for (ip1, ip2, iR, v) in enumerate(int)
		@printf(io, "%6u %6u %6u %14.8f %14.8f\n", ip1, ip2, iR, reim(v)...)
	end

	return nothing
end
