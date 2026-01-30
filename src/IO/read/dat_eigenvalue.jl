
function Base.read(io::IO, ::Type{dat_eigenvalue})

	dat = readdlm(io)
	nband = Int(dat[1])
	dat = dat[2:end]

	nk = Int(length(dat) / nband)
	band = reshape(dat, nband, nk)

	return band
end
