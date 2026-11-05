function Base.read(io::IO, ::Type{dat_eigenvalue})

	data = readdlm(io)
	nband = Int(data[1])
	data = data[2:end]

	nk = Int(length(data) / nband)
	band = reshape(data, nband, nk)

	return band
end
