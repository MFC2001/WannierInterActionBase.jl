"""
	wannierize(band_opt::AbstractVector{<:Eigen}, u_matrix::Array{<:Number, 3}) -> Array{ComplexF64, 3}

Rotate the eigenvectors of every k-point with `u_matrix[:, :, ik]` and return the Wannier
functions as an array of size `(nbasis, nk, nw)`.
"""
function wannierize(band_opt::AbstractVector{<:Eigen}, u_matrix::Array{<:Number, 3})

	(nbasis, nw) = size(band_opt[1].vectors)
	nk = length(band_opt)

	# T = promote_type(eltype(band_opt).parameters[1],eltype(u_matrix))
	wannier = Array{ComplexF64}(undef, nbasis, nk, nw)

	Threads.@threads for k in 1:nk
		for w in 1:nw
			wannier[:, k, w] = band_opt[k].vectors * u_matrix[:, w, k]
		end
	end

	return wannier
end
