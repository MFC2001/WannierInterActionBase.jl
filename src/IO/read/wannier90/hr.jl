
function Base.read(io::IO, ::Type{wannier90_hr}; heps::Real = 0, readimag::Bool = true, μ::Real = 0, hrsort::Bool = false)

	if heps < 0
		error("Please input a positive heps!")
	end

	#first line is comment.
	readline(io)

	norb = parse(Int, readline(io))

	for _ in 1:((parse(Int, readline(io))-1)÷15+1)
		readline(io)
	end
	hr = readdlm(io)

	allpath = Int.(hr[:, 1:5])
	allvalue = hr[:, 6:7]

	if !readimag || all(iszero, allvalue[:, 2])
		allvalue = allvalue[:, 1]
		I = map(x -> abs(x) >= heps, allvalue)

		path = allpath[I, :]
		value = allvalue[I]
	else
		allvalue = complex.(allvalue[:, 1], allvalue[:, 2])
		heps2 = heps^2
		I = map(x -> abs2(x) >= heps2, allvalue)

		path = allpath[I, :]
		value = allvalue[I]
	end

	return wannier90_hr(path, value; orbindex = collect(1:norb), μ, hrsort)
end
