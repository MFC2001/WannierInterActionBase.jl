function Base.read(io::IO, ::Type{vcoul_bgw})
	qpoints = Vector{ReducedCoordinates{Float64}}(undef, 0)
	Glists = Vector{Vector{ReducedCoordinates{Int}}}(undef, 0)
	vals = Vector{Vector{Float64}}(undef, 0)

	data = split(readline(io))
	push!(qpoints, ReducedCoordinates(parse.(Float64, data[1:3])))
	Glist_q = Vector{ReducedCoordinates{Int}}(undef, 0)
	push!(Glist_q, ReducedCoordinates(parse.(Int, data[4:6])))
	push!(Glists, Glist_q)
	vals_q = Vector{Float64}(undef, 0)
	push!(vals_q, parse(Float64, data[7]))
	push!(vals, vals_q)

	for line in readlines(io)
		data = split(line)
		qpoint = ReducedCoordinates(parse.(Float64, data[1:3]))
		if any(≥(1e-6), qpoint - qpoints[end])
			push!(qpoints, qpoint)
			push!(Glists, Vector{ReducedCoordinates{Int}}(undef, 0))
			push!(vals, Vector{Float64}(undef, 0))
		end
		push!(Glists[end], ReducedCoordinates(parse.(Int, data[4:6])))
		push!(vals[end], parse(Float64, data[7]))
	end

	return vcoul_bgw(qpoints, Glists, vals)
end
