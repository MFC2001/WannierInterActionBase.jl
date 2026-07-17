function Base.findfirst(f::Function, z::Iterators.Zip)
	for (idx, val) in enumerate(z)
		f(val) && return idx
	end
	return nothing
end
