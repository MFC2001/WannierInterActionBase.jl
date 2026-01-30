
function Base.write(io::IO, cell::Cell, ::Val{:QEcell}; comment = "")
	now = Dates.format(Dates.now(), "yyyy/m/d H:M:S")
	comment = comment * " " * now

	println(io, comment)
	println(io, "")

	println(io, "CELL_PARAMETERS (angstrom)")
	for i in 1:3
		@printf(io, "%26.16f %26.16f %26.16f\n", cell.lattice[:, i]...)
	end

	println(io, "ATOMIC_POSITIONS (crystal)")
	for i in eachindex(cell.name)
		@printf(io, "%-6s %20.10f %20.10f %20.10f\n", cell.name[i], cell.location[i]...)
	end

	return nothing
end
