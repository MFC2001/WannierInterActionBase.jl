function disentangled_band(band::AbstractVector{<:Eigen}, u_matrix_opt::Array{<:Number, 3}, exclude_bands::AbstractVector{<:Integer})

	if isempty(exclude_bands)
		return disentangled_band(band, u_matrix_opt)
	end

	nband = length(band[1].values)
	bandindex = setdiff(collect(1:nband), exclude_bands)

	band_opt = similar(band)

	if isempty(u_matrix_opt)
		Threads.@threads for k in eachindex(band)
			values = band[k].values[bandindex]
			vectors = band[k].vectors[:, bandindex]
			band_opt[k] = Eigen(values, vectors)
		end
	else
		nwann = size(u_matrix_opt, 2)
		nband = size(u_matrix_opt, 1)
		Threads.@threads for k in eachindex(band)
			values = Vector{Float64}(undef, nwann)
			for i in 1:nwann
				values[i] = sum(n -> abs2(u_matrix_opt[n, i, k]) * band[k].values[bandindex[n]], Base.OneTo(nband))
			end
			vectors = band[k].vectors[:, bandindex] * u_matrix_opt[:, :, k]
			band_opt[k] = Eigen(values, vectors)
		end
	end

	return band_opt
end

function disentangled_band(band::AbstractVector{<:Eigen}, u_matrix_opt::Array{<:Number, 3})

	if isempty(u_matrix_opt)
		band_opt = deepcopy(band)
	else
		band_opt = similar(band)

		nwann = size(u_matrix_opt, 2)
		nband = size(u_matrix_opt, 1)
		Threads.@threads for k in eachindex(band)
			values = Vector{Float64}(undef, nwann)
			for i in 1:nwann
				values[i] = sum(n -> abs2(u_matrix_opt[n, i, k]) * band[k].values[n], Base.OneTo(nband))
			end
			vectors = band[k].vectors * u_matrix_opt[:, :, k]
			band_opt[k] = Eigen(values, vectors)
		end
	end

	return band_opt
end
