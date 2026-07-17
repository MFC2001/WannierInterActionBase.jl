function longrange_ϵ(U_screened::HR, U_bare::HR, lattice, orblocat::AbstractVector{<:AbstractVector{<:Real}}; atol_range = 1e-3, atol_var = 1e-3)
	nU = length(U_bare)
	nU == length(U_screened) || error("Mismatched U!")

	ldot = parent(lattice)
	ldot = Mat3(transpose(ldot) * ldot)
	r_norm = Vector{Float64}(undef, nU)
	for i in Base.OneTo(nU)
		path = view(U_bare.path, :, i)
		R = ReducedCoordinates{Int}(path[1], path[2], path[3])
		r_frac = R + orblocat[path[5]] - orblocat[path[4]]
		r_norm[i] = sqrt(r_frac ⋅ (ldot * r_frac))
	end
	ϵr = real(U_bare.value ./ U_screened.value)

	result = tail_mean(r_norm, ϵr; step = 0.1, rtol_range = 0, atol_range, rtol_var = 0, atol_var)

	return (ϵ = result.mean, rcut = result.xmin, nU = result.num, ϵvar = result.var)
end
function longrange_ϵ(U_screened::HR, U_bare::HR, lattice, orblocat::AbstractVector{<:CartesianCoordinates}; atol_range = 1e-3, atol_var = 1e-3)
	nU = length(U_bare)
	nU == length(U_screened) || error("Mismatched U!")

	lattice = Lattice(parent(lattice))
	orblocat = map(orblocat) do x
		lattice \ x
	end
	return longrange_ϵ(U_screened, U_bare, lattice, orblocat; atol_range, atol_var)
end
function tail_mean(x::AbstractVector{<:Real}, y::AbstractVector{<:Real};
	step::Union{Real, Nothing} = nothing, nstep::Integer = 100,
	rtol_range = 0.1, atol_range = 0,
	rtol_var = 0.2, atol_var = 0,
	minpoints = 5,
)
	n = length(x)
	n == length(y) || error("Mismatched length of x and y!")

	I = sortperm(x)
	x = x[I]
	y = y[I]

	(ymin, ymax) = extrema(y)
	global_range = ymax - ymin
	if global_range < atol_range
		return sum(y)
	end

	global_var = population_variance(y)
	if global_range < atol_var
		return sum(y)
	end

	tol_range = max(rtol_range * global_range, atol_range)
	tol_var = max(rtol_var * global_var, atol_var)

	if isnothing(step)
		nstep > 0 || error("nstep should be positive!")
		step = (x[end] - x[1]) / nstep
	else
		nstep = floor(Int, (x[end] - x[1]) / step)
		nstep > 0 || error("nstep should be positive!")
	end

	left = 0
	n_remain = 0
	final_var = 0.0
	final_range = 0.0
	for istep in 1:nstep
		left = findfirst(≥(x[1] + istep * step), x)
		n_remain = n - left + 1
		if n_remain ≤ minpoints
			@warn "use final $minpoints points!"
			y_remain = view(y, (n-minpoints+1):n)
			return (mean = sum(y_remain) / minpoints, xmin = x[n-minpoints+1], num = minpoints, var = population_variance(y_remain))
		end

		y_remain = view(y, left:n)
		(ymin, ymax) = extrema(y_remain)
		final_range = ymax - ymin
		if final_range ≤ tol_range
			final_var = population_variance(y_remain)
			if final_var ≤ tol_var
				return (mean = sum(y_remain) / n_remain, xmin = x[left], num = n_remain, var = final_var)
			end
		end
	end

	@warn "use final interval!"
	y_remain = view(y, left:n)
	return (mean = sum(y_remain) / n_remain, xmin = x[left], num = n_remain, var = population_variance(y_remain))
end

function population_variance(arr)
	N = length(arr)
	μ = sum(arr) / N
	squared_diff = sum((x - μ)^2 for x in arr)
	return squared_diff / N
end
