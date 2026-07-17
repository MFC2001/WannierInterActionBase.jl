function read_record(io::IO, ::Type{Bool}, dims)
	data = read_record(io)
	data = reinterpret(Int32, data)
	if length(data) ≠ prod(dims)
		error("record length mismatch: expected $(prod(dims)), got $(length(data))")
	end
	return reshape(data, dims) .≠ 0
end
function read_record(io::IO, T::Type, dims)
	data = read_record(io)
	data = reinterpret(T, data)
	if length(data) ≠ prod(dims)
		error("record length mismatch: expected $(prod(dims)), got $(length(data))")
	end
	return reshape(data, dims)
end
function read_record(io::IO, ::Type{String})::String
	data = read_record(io)
	return String(data)
end
function read_record(io::IO, ::Type{Bool})::Bool
	data = read_record(io)
	len_data = length(data)
	len_T = sizeof(Int32)
	if len_data ≠ len_T
		error("record type mismatch: expected $len_T bytes of Int32 for Bool, got $len_data bytes.")
	end
	data = reinterpret(Int32, data)[1]
	return data ≠ 0
end
function read_record(io::IO, T::Type)::T
	data = read_record(io)
	len_data = length(data)
	len_T = sizeof(T)
	if len_data ≠ len_T
		error("record type mismatch: expected $len_T bytes of $T, got $len_data bytes.")
	end
	data = reinterpret(T, data)[1]
	return data
end
function read_record(io::IO, nb::Integer; all = true)::Vector{UInt8}
	data = read_record(io)
	len_data = length(data)
	if len_data ≠ nb
		error("record length mismatch: expected $nb, got $len_data")
	end
	return data
end
function read_record(io::IO)::Vector{UInt8}
	nb_start = read(io, Int32)
	data = read(io, nb_start)
	nb_end = read(io, Int32)
	if nb_start ≠ nb_end
		error("record length mismatch: $nb_start vs $nb_end")
	end
	return data
end
