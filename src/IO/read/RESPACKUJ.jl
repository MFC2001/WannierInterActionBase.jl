
function Base.read(io::IO, ::Type{RESPACKU}; norb::Integer)

	norb2 = norb^2

	#first three lines is comment.
	for _ in 1:3
		readline(io)
	end

	ppath = Vector{Vector{Int}}(undef, 0)
	vvalue = Vector{Matrix{Float64}}(undef, 0)
	while !eof(io)
		mark(io)
		line = readline(io)
		if occursin(r"^\s*(#|$)", line)
			continue
		end

		push!(ppath, parse.(Int, split(line)))

		M = Matrix{Float64}(undef, norb2, 4)
		for i in 1:norb2
			M[i, :] = parse.(Float64, split(readline(io)))
		end
		push!(vvalue, M)
	end

	n = length(ppath)
	path = Matrix{Int}(undef, n * norb2, 5)
	value = Vector{ComplexF64}(undef, n * norb2)
	@views for i in 1:n, (ii, j) in enumerate(norb2*(i-1)+1:norb2*i)
		path[j, :] = [ppath[i]; vvalue[i][ii, 1:2]]
		value[j] = complex.(vvalue[i][ii, 3], vvalue[i][ii, 4])
	end

	return HR(path, value; orbindex = collect(1:norb), hrsort = 'N')
end

function Base.read(io::IO, ::Type{RESPACKJ}; norb::Integer)

	norb2 = norb^2

	#first three lines is comment.
	for _ in 1:3
		readline(io)
	end

	ppath = Vector{Vector{Int}}(undef, 0)
	vvalue = Vector{Matrix{Float64}}(undef, 0)
	while !eof(io)
		mark(io)
		line = readline(io)
		if occursin(r"^\s*(#|$)", line)
			continue
		end

		push!(ppath, parse.(Int, split(line)))

		M = Matrix{Float64}(undef, norb2, 4)
		for i in 1:norb2
			M[i, :] = parse.(Float64, split(readline(io)))
		end
		push!(vvalue, M)
	end

	n = length(ppath)
	path = Matrix{Int}(undef, n * norb2, 5)
	value = Vector{ComplexF64}(undef, n * norb2)
	@views for i in 1:n
		if iszero(ppath[i])
			for (ii, j) in enumerate(norb2*(i-1)+1:norb2*i)
				path[j, :] .= [ppath[i]; vvalue[i][ii, 1:2]]
				if path[j, 4] == path[j, 5]
					value[j] = 0
				else
					value[j] = complex.(vvalue[i][ii, 3], vvalue[i][ii, 4])
				end
			end
		else
			for (ii, j) in enumerate(norb2*(i-1)+1:norb2*i)
				path[j, :] = [ppath[i]; vvalue[i][ii, 1:2]]
				value[j] = complex.(vvalue[i][ii, 3], vvalue[i][ii, 4])
			end
		end
	end

	return HR(path, value; orbindex = collect(1:norb), hrsort = 'N')
end
