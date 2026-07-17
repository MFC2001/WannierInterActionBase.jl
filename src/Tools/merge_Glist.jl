function merge_Glist(Glists::Vector{ReducedCoordinates{Int}}...;
	rldot::Union{AbstractMatrix{<:Real}, Nothing} = nothing, Gidx::Bool = true)

	#create Glist_unique
	#flatten will iterate every elements in a nested structure.
	Glist_unique = Set(Iterators.flatten(Glists))
	Glist_unique = collect(Glist_unique)

	#sort by norm
	if isnothing(rldot)
		sort_keys = map(Glist_unique) do G
			G ⋅ G
		end
	else
		sort_keys = map(Glist_unique) do G
			G ⋅ (rldot * G)
		end
	end
	perm = sortperm(sort_keys)
	permute!(Glist_unique, perm)
	# sort!(Glist_unique, by = G -> sum(abs2, G))

	if Gidx
		G_to_idx = Dict{ReducedCoordinates{Int}, Int}(G => i for (i, G) in enumerate(Glist_unique))
		Gidx = Vector{Vector{Int}}(undef, length(Glists))
		for i in eachindex(Glists)
			Gidx[i] = map(G -> G_to_idx[G], Glists[i])
		end
		return Glist_unique, Gidx
	else
		return Glist_unique
	end
end
