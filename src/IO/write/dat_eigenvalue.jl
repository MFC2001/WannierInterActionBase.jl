
function Base.write(io::IO, band::Matrix{<:Real}, ::Type{dat_eigenvalue})

	println(io, size(band, 1))
	writedlm(io, band[:])

	return nothing
end
