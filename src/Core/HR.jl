export HR, numorb, prune, shift_fermi_energy!, spinHR, translate, translate!, reindexHR, reindexHR!
"""
	HR{T <: Number}

The data in HR.dat. Fields:
- `orbindex::Vector{Int}`: the index of orbtal, usually equal to 1:norb;
- `path::Matrix{Int}`: its dimension is (5, N), each column is [Rx, Ry, Rz, i, j];
- `value::Vector{T}`: its length is N.
"""
struct HR{T <: Number}
	orbindex::Vector{Int}
	path::Matrix{Int} # [5, N], one row is one hopping.
	value::Vector{T}
end
"""
	HR(path, value::Vector{T}; orbindex = sort(unique(path[4:5, :])), μ::Real = 0, hrsort::Bool = false) -> HR{T}

Construct a `HR` object from the given hopping path indices and corresponding numerical values.

# Arguments
- `path`: An integer matrix with the size of (5, N), storing the hopping path indices.
- `value`: A vector of length N, storing the numerical values corresponding to the hopping paths in `path`.
- `orbindex`: Optional; A list of orbital indices, with the default value (`sort(unique(path[4:5, :]))`).
- `μ`: Optional; The Fermi energy shift parameter, with a default value of 0.
- `hrsort`: Optional; A flag to control whether to sort the hopping data by the `path` indices.
"""
function HR(path::AbstractMatrix{<:Integer}, value::AbstractVector{T};
	orbindex::AbstractVector{<:Integer} = sort(unique(path[4:5, :])), μ::Real = 0, hrsort::Bool = false)::HR{T} where {T <: Number}

	path_rows, path_cols = size(path)
	value_len = length(value)
	@assert path_rows == 5 "Invalid `path` matrix dimensions: expected 5 rows, got $path_rows rows."
	@assert path_cols == value_len "Dimension mismatch between `path` and `value`: `path` has $path_cols columns, but `value` has $value_len elements."

	hr = HR{T}(orbindex, path, value)

	if abs(μ) > 1e-6
		shift_fermi_energy!(hr, μ)
	end

	if hrsort[1] ∈ ['Y', 'y']
		sort!(hr)
	end

	return hr
end

function Base.show(io::IO, hr::HR)
	print(io, "HR with $(length(hr)) hoppings and $(numorb(hr)) orbitals.")
end
Base.length(hr::HR) = length(hr.value)
numorb(hr::HR) = length(hr.orbindex)
Base.convert(::Type{HR{T₁}}, hr::HR{T₂}) where {T₁, T₂} =
	HR{T₁}(hr.orbindex, hr.path, convert(Vector{T₁}, hr.value))

"""
	filter(F, hr::HR{T}) where {T}

F(i::Int, j::Int, R::Vector{Int}, value::T) -> Bool.
t_{i0,jR} = value.
"""
function Base.filter(F, hr::HR{T}) where {T}

	I = similar(hr.value, Bool)
	for i in eachindex(I)
		I[i] = F(hr.path[4, i], hr.path[5, i], hr.path[1:3, i], hr.value[i])
	end

	new_path = hr.path[:, I]
	new_value = hr.value[I]

	return HR{T}(hr.orbindex, new_path, new_value)
end
function prune(hr::HR{T}, cutoff::Real) where {T <: Real}
	@assert cutoff ≥ 0 "cutoff should be positive!"
	f(ip1, ip2, iR, v) = abs(v) ≥ cutoff
	return filter(f, hr)
end
function prune(hr::HR{T}, cutoff::Real) where {T <: Complex}
	@assert cutoff ≥ 0 "cutoff should be positive!"
	cutoff2 = cutoff^2
	f(ip1, ip2, iR, v) = abs2(v) ≥ cutoff2
	return filter(f, hr)
end
"""
	sort(hr::HR{T}) where {T}

Sort hr.path and hr.value by hr.path.
"""
function Base.sort(hr::HR{T}) where {T}
	pathI = [hr.path; reshape(collect(eachindex(hr.value)), 1, :)]
	pathI = sortslices(pathI; dims = 2)
	I = pathI[6, :]

	sort_path = pathI[1:5, :]
	sort_value = hr.value[I]

	return HR{T}(hr.orbindex, sort_path, sort_value)
end
"""
	sort!(hr::HR{T}) where {T}

Sort hr.path and hr.value by hr.path.
"""
function Base.sort!(hr::HR{T}) where {T}
	pathI = [hr.path; reshape(collect(eachindex(hr.value)), 1, :)]
	pathI = sortslices(pathI; dims = 2)
	I = pathI[6, :]

	hr.path .= pathI[1:5, :]

	value_copy = copy(hr.value)
	hr.value .= value_copy[I]

	return hr
end
"""
	shift_fermi_energy!(hr::HR{T}, μ::Real)::HR{T} where{T}
"""
function shift_fermi_energy!(hr::HR{T}, μ::Real)::HR{T} where {T}
	I = findall(x -> iszero(x[1:3]) && isequal(x[4], x[5]), eachcol(hr.path))
	hr.value[I] .-= μ
	return hr
end
function Base.union(hr₁::HR{T₁}, hr₂::HR{T₂}) where {T₁, T₂}
	Set(hr₁.orbindex) == Set(hr₂.orbindex) || error("Mismatched hr.")
	path = [hr₁.path hr₂.path]
	value = [hr₁.value; hr₂.value]
	return unique(HR{promote_type(T₁, T₂)}(hr₁.orbindex, path, value))
end
function Base.unique(hr::HR{T})::HR{T} where {T}
	unique_path2value = Dict{SVector{5, Int}, T}()
	for (path, value) in zip(eachcol(hr.path), hr.value)
		path = SVector{5, Int}(path)
		if haskey(unique_path2value, path)
			unique_path2value[path] += value
		else
			unique_path2value[path] = value
		end
	end
	n = length(unique_path2value)
	unique_path = Matrix{Int}(undef, 5, n)
	unique_value = Vector{T}(undef, n)
	i = 0
	for (key, value) in unique_path2value
		i += 1
		unique_path[:, i] = key
		unique_value[i] = value
	end
	unique_hr = HR{T}(hr.orbindex, unique_path, unique_value)
	return sort!(unique_hr)
end

"""
	spinHR(hr::HR{T}; mode = conj)::HR{T} where {T}

Generate a spinful hr from a spinless hr.
Acctually, `spinHR` just makes a copy, it won't introduce new hopping.
"""
function spinHR(hr::HR{T}; mode = conj)::HR{T} where {T}
	up_path = hr.path
	up_value = hr.value

	norb = length(hr.orbindex)
	#make sure using a seperate path
	dn_path = copy(hr.path)
	dn_path[:, 4:5] .+= norb
	#conj is fron time-reversal.
	dn_value = mode.(hr.value)

	sum_path = [up_path; dn_path]
	sum_value = [up_value; dn_value]

	sum_orbindex = sort([hr.orbindex; hr.orbindex .+ norb])

	return HR{T}(sum_orbindex, sum_path, sum_value)
end

"""
	translate(hr::HR{T}, centres::Pair{<:Integer, <:AbstractVector{<:Integer}}...) -> HR{T}

You only need to input the wannier centres R which isn't in the [0, 0, 0] cell.

See also `translate!`.
"""
function translate(hr::HR{T}, aim_centres::Pair{<:Integer, <:AbstractVector{<:Integer}}...)::HR{T} where {T}
	nw = maximum(hr.orbindex)
	centres = [ReducedCoordinates{Int}(0, 0, 0) for _ in 1:nw]
	for (iw, R) in aim_centres
		centres[iw] = ReducedCoordinates{Int}(R)
	end

	hr_path_new = similar(hr.path)
	for ip in axes(hr.path, 2)
		R = ReducedCoordinates{Int}(hr.path[1:3, ip])
		iw1, iw2 = hr.path[4, ip], hr.path[5, ip]
		hr_path_new[4, ip] = iw1
		hr_path_new[5, ip] = iw2
		hr_path_new[1:3, ip] .= (R + centres[iw2] - centres[iw1])
	end
	return HR{T}(copy(hr.orbindex), hr_path_new, copy(hr.value))
end
"""
	translate!(hr::HR, centres::Pair{<:Integer, <:AbstractVector{<:Integer}}...)

You only need to input the wannier centres R which isn't in the [0, 0, 0] cell.

See also `translate`.
"""
function translate!(hr::HR{T}, aim_centres::Pair{<:Integer, <:AbstractVector{<:Integer}}...)::HR{T} where {T}
	nw = maximum(hr.orbindex)
	centres = [ReducedCoordinates{Int}(0, 0, 0) for _ in 1:nw]
	for (iw, R) in aim_centres
		centres[iw] = ReducedCoordinates{Int}(R)
	end
	for ip in axes(hr.path, 2)
		R = ReducedCoordinates{Int}(hr.path[1:3, ip])
		iw1, iw2 = hr.path[4, ip], hr.path[5, ip]
		hr.path[1:3, ip] .= (R + centres[iw2] - centres[iw1])
	end
	return hr
end

"""
	reindexHR(hr::HR{T}, newindex::AbstractVector{<:Integer})::HR{T}

This function will change the orbindex and orbpath of hr, and don't need the index of hr.
"""
function reindexHR(hr::HR{T}, newindex::AbstractVector{<:Integer})::HR{T} where {T}

	if numorb(hr) ≠ length(newindex)
		error("Wrong length of newindex to hr.")
	end

	neworbpath = Matrix{Int}(undef, 2, length(hr))
	@views for i in eachindex(hr.orbindex)
		I = hr.path[4:5, :] .== hr.orbindex[i]
		neworbpath[I] .= newindex[i]
	end

	newpath = [hr.path[1:3, :]; neworbpath]

	return HR{T}(newindex, newpath, copy(hr.value))
end
"""
	reindexHR(hr::HR, N::Integer)::HR
	reindexHR!(hr::HR, N::Integer)::HR

This function will change the orbindex and orbpath of hr, and add N to these two data.
"""
function reindexHR(hr::HR{T}, N::Integer)::HR{T} where {T}
	return reindexHR!(deepcopy(hr), N)
end
function reindexHR!(hr::HR{T}, N::Integer)::HR{T} where {T}
	hr.path[4:5, :] .+= N
	return hr
end
