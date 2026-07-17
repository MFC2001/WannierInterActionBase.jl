function Base.write(io::IO, hr::HR, ::Type{wannier90_hr}; comment = "")
	now = Dates.format(Dates.now(), "yyyy/m/d H:M:S")
	comment = comment * ", " * now

	println(io, comment)

	@printf(io, "%12u\n", numorb(hr))
	allpath = unique(x[1:3] for x in eachcol(hr.path))
	nkpoint = length(allpath)
	@printf(io, "%12u\n", nkpoint)

	for i in 1:nkpoint
		@printf(io, "%5u", 1)
		if mod(i, 15) == 0
			@printf(io, "\n")
		end
	end
	if mod(nkpoint, 15) ≠ 0
		@printf(io, "\n")
	end

	for i in eachindex(hr.value)
		@printf(io, "%5u %5u %5u %7u %7u %14.8f %14.8f\n", hr.path[:, i]..., reim(hr.value[i])...)
	end

	return nothing
end
