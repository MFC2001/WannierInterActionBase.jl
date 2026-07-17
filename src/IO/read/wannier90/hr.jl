function Base.read(io::IO, ::Type{wannier90_hr};
	atol::Real = 0, readimag::Bool = false, μ::Real = 0, hrsort::Bool = false)

	@assert atol ≥ 0 "atol should be positive!"

	#first line is comment.
	readline(io)

	norb = parse(Int, readline(io))

	nR = parse(Int, readline(io))
	nline = (nR - 1) ÷ 15
	ndegen = Vector{Int}(undef, nR)
	for i in 1:nline
		ndegen[((i-1)*15+1):(i*15)] .= parse.(Int, split(readline(io)))
	end
	ndegen[(nline*15+1):end] .= parse.(Int, split(readline(io)))

	hr = readdlm(io)

	allpath = transpose(Int.(hr[:, 1:5]))
	allvalue = hr[:, 6:7]

	if !readimag && all(iszero, allvalue[:, 2])
		allvalue = allvalue[:, 1]
		I = map(x -> abs(x) ≥ atol, allvalue)

		path = allpath[:, I]
		value = allvalue[I]
	else
		allvalue = complex.(allvalue[:, 1], allvalue[:, 2])
		atol2 = atol^2
		I = map(x -> abs2(x) ≥ atol2, allvalue)

		path = allpath[:, I]
		value = allvalue[I]
	end

	npath = size(path, 2)
	Rlist = Vector{ReducedCoordinates{Int}}(undef, npath)
	Rlist_unique = Vector{ReducedCoordinates{Int}}(undef, 0)
	for i in 1:npath
		R = ReducedCoordinates{Int}(path[1, i], path[2, i], path[3, i])
		Rlist[i] = R
		if R ∉ Rlist_unique
			push!(Rlist_unique, R)
		end
	end
	if length(Rlist_unique) ≠ nR
		error("The `nR` in hr file is wrong!")
	end

	R_to_idx = Dict([R => i for (i, R) in enumerate(Rlist_unique)])

	for i in 1:npath
		idx = R_to_idx[Rlist[i]]
		value[i] /= ndegen[idx]
	end

	return HR(path, value; orbindex = collect(1:norb), μ, hrsort)
end
