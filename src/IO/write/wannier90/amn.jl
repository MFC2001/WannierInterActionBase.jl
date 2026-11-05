function Base.write(io::IO, A::AbstractArray{<:Number}, ::Type{wannier90_amn}; comment = "")
	now = Dates.format(Dates.now(), "yyyy/m/d H:M:S")
	comment = comment * " " * now

	(nband, nw, nk) = size(A)

	println(io, comment)

	@printf(io, "%12u %12u %12u \n", nband, nk, nw)

	for k in 1:nk
		for iw in 1:nw
			for n in 1:nband
				@printf(io, "%5u %5u %5u %18.12f %18.12f \n", n, iw, k, reim(A[n, iw, k])...)
			end
		end
	end

	return nothing
end
