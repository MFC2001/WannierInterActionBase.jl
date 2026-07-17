
function Base.read(io::IO, ::Type{WanIntijkl0}; readimag::Bool = true)

	#first line is comment.
	readline(io)

	nw = parse(Int, readline(io))
	data = readdlm(io)
	@assert size(data, 1) == nw^4 "Wrong file of interaction between wannier function in 0-cell."

	data_idx = Int.(data[:, 1:4])
	if !readimag && all(iszero, intvalue[:, 2])
		data_value = data[:, 5]
		intvalue = Array{Float64}(undef, nw, nw, nw, nw)
	else
		data_value = complex.(data[:, 5], data[:, 6])
		intvalue = Array{ComplexF64}(undef, nw, nw, nw, nw)
	end

	for idx in 1:(nw^4)
		i, j, k, l = data_idx[idx, :]
		intvalue[i, j, k, l] = data_value[idx]
	end

	return intvalue
end

function Base.read(io::IO, ::Type{WanIntijklR}; readimag::Bool = false, atol::Real = 1e-6)
	#first line is comment.
	readline(io)

	(nw, np, nR, nv) = parse.(Int, split(readline(io)))

	# pair
	readline(io)
	pair = Vector{NTuple{5, Int}}(undef, np)
	for ip in 1:np
		(_, iw1, iw2, R1, R2, R3) = parse.(Int, split(readline(io)))
		pair[ip] = (iw1, iw2, R1, R2, R3)
	end

	pR_all = map(pair) do p
		ReducedCoordinates{Int}(p[3], p[4], p[5])
	end
	pR = unique(pR_all)
	npR = length(pR)
	pR_to_idx = Dict{ReducedCoordinates{Int}, Int}(R => i for (i, R) in enumerate(pR))
	pair = map(pair) do p
		(p[1], p[2], pR_to_idx[ReducedCoordinates{Int}(p[3], p[4], p[5])])
	end

	# R
	readline(io)
	R = Vector{ReducedCoordinates{Int}}(undef, nR)
	for iR in 1:nR
		(_, R1, R2, R3) = parse.(Int, split(readline(io)))
		R[iR] = ReducedCoordinates{Int}(R1, R2, R3)
	end

	#value
	readline(io)
	ips = Vector{Int}(undef, nv)
	jps = Vector{Int}(undef, nv)
	iRs = Vector{Int}(undef, nv)
	value = Vector{ComplexF64}(undef, nv)
	allisreal = true
	for i in Base.OneTo(nv)
		line = split(readline(io))
		ips[i] = parse(Int, line[1])
		jps[i] = parse(Int, line[2])
		iRs[i] = parse(Int, line[3])
		real_part = parse(Float64, line[4])
		imag_part = parse(Float64, line[5])
		value[i] = complex(real_part, imag_part)
		allisreal &= (abs(imag_part) < atol)
	end

	if !readimag && allisreal
		value = real(value)
	end

	int = WanIntijklR(nw, npR, pR, np, pair, nR, R, nv, ips, jps, iRs, value)

	return prune(int, atol)
end
