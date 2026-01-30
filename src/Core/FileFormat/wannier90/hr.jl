export wannier90_hr, numorb, shift_fermi_energy!, spinhr, reindexhr, reindexhr!
"""
wannier90_hr.dat

	read(path/io, ::Type{wannier90_hr}; heps::Real = 0, readimag::Bool = true, μ::Real = 0, hrsort::Bool = false) -> wannier90_hr

- `heps`: can be set as a random positive real number, which is used to filter small values.
- `readimag`: decides whether to read imaginary part.
- `μ`: can be set as a random real number, which is the chemical potential.
- `hrsort`: decides whether to sort the object according to the standard order.

	write(path/io, hr::wannier90_hr, ::Type{wannier90_hr})

Write hr to `path/io` in the format of `wannier90_hr.dat`.

	wannier90_hr{T <: Number} <: FileFormat

The data in wannier90_hr.dat. Fields:
- `orbindex::Vector{Int}`: the index of orbtal, usually equal to 1:norb;
- `path::Matrix{Int}`: its dimension is (5, N), each column is [Rx, Ry, Rz, i, j];
- `value::Vector{T}`: its length is N.

"""
struct wannier90_hr{T <: Number} <: FileFormat
	orbindex::Vector{Int}
	path::Matrix{Int} # [5, N], one row is one hopping.
	value::Vector{T}
end
"""
	wannier90_hr(path, value::Vector{T}; orbindex = sort(unique(path[4:5, :])), μ::Real = 0, hrsort::Bool = false) -> wannier90_hr{T}

Construct a `wannier90_hr` object from the given hopping path indices and corresponding numerical values.

# Arguments
- `path`: An integer matrix with the size of (5, N), storing the hopping path indices.
- `value`: A vector of length N, storing the numerical values corresponding to the hopping paths in `path`.
- `orbindex`: Optional; A list of orbital indices, with the default value (`sort(unique(path[4:5, :]))`).
- `μ`: Optional; The Fermi energy shift parameter, with a default value of 0.
- `hrsort`: Optional; A flag to control whether to sort the hopping data by the `path` indices.
"""
function wannier90_hr(path::AbstractMatrix{<:Integer}, value::AbstractVector{T};
	orbindex::AbstractVector{<:Integer} = sort(unique(path[4:5, :])), μ::Real = 0, hrsort::Bool = false)::wannier90_hr{T} where {T <: Number}

	path_rows, path_cols = size(path)
	value_len = length(value)
	@assert path_rows == 5 "Invalid `path` matrix dimensions: expected 5 rows, got $path_rows rows."
	@assert path_cols == value_len "Dimension mismatch between `path` and `value`: `path` has $path_cols columns, but `value` has $value_len elements."

	hr = wannier90_hr{T}(orbindex, path, value)

	if abs(μ) > 1e-6
		shift_fermi_energy!(hr, μ)
	end

	if hrsort[1] ∈ ['Y', 'y']
		sort!(hr)
	end

	return hr
end

function Base.show(io::IO, hr::wannier90_hr)
	print(io, "wannier90_hr with $(length(hr)) hoppings and $(numorb(hr)) orbitals.")
end
Base.length(hr::wannier90_hr) = length(hr.value)
numorb(hr::wannier90_hr) = length(hr.orbindex)

Base.convert(::wannier90_hr{T₁}, hr::wannier90_hr{T₂}) where {T₁, T₂} =
	wannier90_hr{T₁}(hr.orbindex, hr.path, convert(Vector{T₁}, hr.value))
Base.convert(::Type{wannier90_hr{T₁}}, hr::wannier90_hr{T₂}) where {T₁, T₂} =
	wannier90_hr{T₁}(hr.orbindex, hr.path, convert(Vector{T₁}, hr.value))

"""
	filter(F, hr::wannier90_hr{T}) where {T}

F(i, j, R, value) -> Bool.
t_{i0,jR} = value.
"""
function Base.filter(F, hr::wannier90_hr{T}) where {T}

	I = similar(hr.value, Bool)
	for i in eachindex(I)
		I[i] = F(hr.path[4, i], hr.path[5, i], hr.path[1:3, i], hr.value[i])
	end

	new_path = hr.path[:, I]
	new_value = hr.value[I]

	return wannier90_hr{T}(hr.orbindex, new_path, new_value)
end
"""
	sort(hr::wannier90_hr{T}) where {T}

Sort hr.path and hr.value by hr.path.
"""
function Base.sort(hr::wannier90_hr{T}) where {T}
	pathI = [hr.path; reshape(collect(eachindex(hr.value)), 1, :)]
	pathI = sortslices(pathI; dims = 2)
	I = pathI[6, :]

	sort_path = pathI[1:5, :]
	sort_value = hr.value[I]

	return wannier90_hr{T}(hr.orbindex, sort_path, sort_value)
end
"""
	sort!(hr::wannier90_hr{T}) where {T}

Sort hr.path and hr.value by hr.path.
"""
function Base.sort!(hr::wannier90_hr{T}) where {T}
	pathI = [hr.path; reshape(collect(eachindex(hr.value)), 1, :)]
	pathI = sortslices(pathI; dims = 2)
	I = pathI[6, :]

	hr.path .= pathI[1:5, :]

	value_copy = copy(hr.value)
	hr.value .= value_copy[I]

	return hr
end
"""
	shift_fermi_energy!(hr::wannier90_hr{T}, μ::Real)::wannier90_hr{T} where{T}
"""
function shift_fermi_energy!(hr::wannier90_hr{T}, μ::Real)::wannier90_hr{T} where {T}
	I = findall(x -> iszero(x[1:3]) && isequal(x[4], x[5]), eachcol(hr.path))
	hr.value[I] .-= μ
	return hr
end
function Base.union(hr₁::wannier90_hr{T₁}, hr₂::wannier90_hr{T₂}) where {T₁, T₂}
	Set(hr₁.orbindex) == Set(hr₂.orbindex) || error("Mismatched hr.")
	path = [hr₁.path hr₂.path]
	value = [hr₁.value; hr₂.value]
	return unique(wannier90_hr{promote_type(T₁, T₂)}(hr₁.orbindex, path, value))
end
function Base.unique(hr::wannier90_hr{T})::wannier90_hr{T} where {T}
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
	unique_hr = wannier90_hr{T}(hr.orbindex, unique_path, unique_value)
	return sort!(unique_hr)
end

"""
	spinhr(hr::wannier90_hr{T}; mode = conj)::wannier90_hr{T} where {T}

Generate a spinful hr from a spinless hr.
Acctually, `spinhr` just makes a copy, it won't introduce new hopping.
"""
function spinhr(hr::wannier90_hr{T}; mode = conj)::wannier90_hr{T} where {T}
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

	return wannier90_hr{T}(sum_orbindex, sum_path, sum_value)
end

"""
	reindexhr(hr::wannier90_hr{T}, newindex::AbstractVector{<:Integer})::wannier90_hr{T}

This function will change the orbindex and orbpath of hr, and don't need the index of hr.
"""
function reindexhr(hr::wannier90_hr{T}, newindex::AbstractVector{<:Integer})::wannier90_hr{T} where {T}

	if numorb(hr) ≠ length(newindex)
		error("Wrong length of newindex to hr.")
	end

	neworbpath = Matrix{Int}(undef, 2, length(hr))
	@views for i in eachindex(hr.orbindex)
		I = hr.path[4:5, :] .== hr.orbindex[i]
		neworbpath[I] .= newindex[i]
	end

	newpath = [hr.path[1:3, :]; neworbpath]

	return wannier90_hr{T}(newindex, newpath, deepcopy(hr.value))
end
"""
	reindexhr(hr::wannier90_hr, N::Integer)::wannier90_hr
	reindexhr!(hr::wannier90_hr, N::Integer)::wannier90_hr

This function will change the orbindex and orbpath of hr, and add N to these two data.
"""
function reindexhr(hr::wannier90_hr{T}, N::Integer)::wannier90_hr{T} where {T}
	return reindexhr!(deepcopy(hr), N)
end
function reindexhr!(hr::wannier90_hr{T}, N::Integer)::wannier90_hr{T} where {T}
	hr.path[4:5, :] .+= N
	return hr
end
