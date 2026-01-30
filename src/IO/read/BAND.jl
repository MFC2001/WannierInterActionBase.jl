
function Base.read(io::IO, ::Type{BAND_dat}; readkname = 'N', smoothband = 'N')
	#first line is comment.
	readline(io)

	line = readline(io)
	line = split(line, ':')[2]
	line = parse.(Int, split(line))
	Nk = line[1]
	Nband = line[2]

	#when output a band.dat by write in this package, there will be a block of kline on the head of band.dat.
	if readkname[1] ∈ ['N', 'n']
		kindex = Int[]
		kname = String[]
	elseif readkname[1] ∈ ['Y', 'y']
		#kindex.
		line = readline(io)
		line = split(line, ':')[2]
		kindex = parse.(Int, split(line))
		#kpoint
		readline(io)
		#kname.
		line = readline(io)
		line = split(line, ':')[2]
		kname = split(line)
		readline(io)
	end

	T = readdlm(io; comments = true, comment_char = '#')

	kline = sort(T[1:Nk, 1])
	band = Matrix{eltype(T)}(undef, Nband, Nk)
	for i in 1:Nband
		index = sortperm(T[(i-1)*Nk+1:i*Nk, 1])
		band[i, index] = T[(i-1)*Nk+1:i*Nk, 2]
	end

	if smoothband[1] ∈ ['Y', 'y']
		smoothband(band)
	end

	return band, Kline(; line = kline, name = kname, index = kindex)
end
"""
	smoothband(band::AbstractMatrix{<:Real})::AbstractMatrix{<:Real}

Just as its name implies.
"""
function smoothband(band::AbstractMatrix{<:Real})::AbstractMatrix{<:Real}

	band = deepcopy(band)
	Nk = size(band, 2)

	for i in 2:Nk
		times = 0 #A counter for avoiding dead loops.
		while true
			times += 1

			dE = band[:, i] - band[:, i-1]
			dE = dE[.!isnan.(dE)]

			(Emax, index) = findmax(abs, dE)
			if Emax < 0.2
				break
			end
			if dE[index] > 0
				band[:, i] = [NaN; band[1:end-1, i]]
			elseif dE[index] < 0
				band[:, i] = [band[2:end, i]; NaN]
			end

			if times > 100
				error("Wrong band from smoothband!")
			end
		end
	end

	I = .!any(ϵ -> isnan(ϵ), band, dims = 2)[:]
	band = band[I, :]

	return band
end



"""
	wannierband(file::AbstractString, nk::Integer)

Read wannier90_band.dat
To be continued.
"""
function wannierband(file::AbstractString, nk::Integer)

	error("To be continued.")
	file = open(file, "r")

	kline = Vector{Float64}(undef, nk)
	E = Matrix{Float64}(undef, 0, nk)

	while !eof(file)
		dE = Vector{Float64}(undef, nk)
		for i in 1:nk
			(kline[i], dE[i]) = parse.(Float64, split(readline(file)))
		end
		E = [E; transpose(dE)]
		readline(file)
	end
	close(file)

	return E, kline
end
