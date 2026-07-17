
function Base.read(io::IO, ::Type{wannier90_UNK}; convert = "native")

	io = FortranFile(io; convert)

	(ngx, ngy, ngz, ik, nbnd) = read(io, (Int32, 5))

	u = Vector{Array{ComplexF64, 3}}(undef, nbnd)
	for i in 1:nbnd
		u[i] = read(io, (ComplexF64, ngx, ngy, ngz))
	end

	return wannier90_UNK(ik, nbnd, (ngx, ngy, ngz), u)
end
