function unique_groups(data)
	d = Dict{eltype(data), Vector{Int}}()
	for (i, x) in enumerate(data)
		push!(get!(d, x, Int[]), i)
	end
	ks = Vector{eltype(data)}()
	vs = Vector{Vector{Int}}()
	sizehint!(ks, length(d))
	sizehint!(vs, length(d))
	for (k, v) in d
		push!(ks, k)
		push!(vs, v)
	end
	return ks, vs
end
function unique_groups(F, data)
	d = Dict{eltype(data), Vector{Int}}()
	for (i, x) in enumerate(data)
		push!(get!(d, F(x), Int[]), i)
	end
	ks = Vector{eltype(data)}()
	vs = Vector{Vector{Int}}()
	sizehint!(ks, length(d))
	sizehint!(vs, length(d))
	for (k, v) in d
		push!(ks, k)
		push!(vs, v)
	end
	return ks, vs
end
