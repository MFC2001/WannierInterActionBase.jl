function Base.show(io::IO, ::MIME"text/plain", lattice::Lattice)
	summary(io, lattice)
	println(io)
	join(io, ' ' * join(row, "  ") * '\n' for row in eachrow(parent(lattice)))
	return nothing
end
function Base.show(io::IO, ::MIME"text/plain", lattice::ReciprocalLattice)
	summary(io, lattice)
	println(io)
	join(io, ' ' * join(row, "  ") * '\n' for row in eachrow(parent(lattice)))
	return nothing
end
# function Base.show(io::IO, ::MIME"text/plain", cell::Cell)
function Base.show(io::IO, cell::Cell)
	summary(io, cell)
	println(io)
	println(io, " lattice:")
	for row in eachrow(parent(Lattice(cell)))
		@printf(io, "%13.5f %13.5f %13.5f\n", row...)
	end
	num_atom = numatom(cell)
	println(io, " $num_atom atoms:")
	if isempty(cell.name)
		for x in cell.location
			@printf(io, "%13.5f %13.5f %13.5f\n", x...)
		end
	else
		for (name, x) in zip(cell.name, cell.location)
			@printf(io, "  %s %13.5f %13.5f %13.5f\n", name, x...)
		end
	end
	return nothing
end