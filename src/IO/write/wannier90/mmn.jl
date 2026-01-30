
function Base.write(io::IO, M::AbstractArray{<:Number}, ::Type{wannier90_mmn}; comment = "",
	nnkpts::AbstractMatrix{<:Integer})

	now = Dates.format(Dates.now(), "yyyy/m/d H:M:S")
	comment = comment * " " * now

	println(io, comment)

	nband = size(M, 1)
	Nk = maximum(nnkpts[1, :])
	nb = Int(size(nnkpts, 2) // Nk)

	@printf(io, "%12u %12u %12u \n", nband, Nk, nb)

	for i in axes(M, 3)
		@printf(io, "%5u %5u %5u %5u %5u \n", nnkpts[:, i]...)
		for n in Base.OneTo(nband)
			for m in Base.OneTo(nband)
				@printf(io, "%18.13f %18.13f \n", reim(M[m, n, i])...)
			end
		end
	end

	return nothing
end
