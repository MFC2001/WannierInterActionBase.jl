
function Base.write(io::IO, hr::wannier90_hr, ::Type{wannier90_hr}; comment = "")
	now = Dates.format(Dates.now(), "yyyy/m/d H:M:S")
	comment = comment * " " * now

	println(io, comment)

	@printf(io, "%12u\n", numorb(hr))
	allpath = union(x[1:3] for x in eachrow(hr.path))
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
		@printf(io, "%5u %5u %5u %7u %7u %12.6f %12.6f\n", hr.path[i, :]..., reim(hr.value[i])...)
	end

	return nothing
end

# Base.write(io::IO, TB::AbstractTightBindModel, ::Type{wannier90_hr}; kwargs...) =
# 	write(io, HR(TB), wannier90_hr; kwargs...)
