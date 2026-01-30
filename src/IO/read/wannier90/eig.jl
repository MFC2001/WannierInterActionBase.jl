
"""
	read(path::AbstractString/io::IO, ::Type{wannier90_eig})

Read wannier90.eig, return a matrix with size (Nband, Nk).
"""
function Base.read(io::IO, ::Type{wannier90_eig})

	data = readdlm(io)

	nband = Int(data[end, 1])
	nk = Int(data[end, 2])

	band = data[:, 3]

	return reshape(band, nband, nk)
end
