function Base.write(io::IO, cell::Cell, ::Type{QEpw}; comment = "",
	title::AbstractString = "",
	kpoints::Union{AbstractBrillouinZone, AbstractVector{<:AbstractVector{<:Real}}, Nothing} = nothing,
)
	now = Dates.format(Dates.now(), "yyyy/m/d H:M:S")
	comment = "! " * comment * " " * now

	println(io, comment)
	println(io, "&CONTROL")
	println(io, "  title = ", title)
	println(io, "/\n")


	println(io, "CELL_PARAMETERS (angstrom)")
	for i in 1:3
		@printf(io, "%26.16f %26.16f %26.16f\n", cell.lattice[:, i]...)
	end
	println(io, "ATOMIC_POSITIONS (crystal)")
	for i in eachindex(cell.name)
		@printf(io, "%-6s %20.10f %20.10f %20.10f\n", cell.name[i], cell.location[i]...)
	end
	println(io, "")

	isnothing(kpoints) || begin
		if kpoints isa MonkhorstPack
			kpoints = RedKgrid(kpoints).kdirect
		elseif kpoints isa RedKgrid
			kpoints = kpoints.kdirect
		elseif kpoints isa Kline
			kpoints = kpoints.kdirect
		elseif kpoints isa AbstractVector{<:AbstractVector{<:Real}}
		else
			error("Wrong grid!")
		end
		println(io, "K_POINTS crystal")
		@printf(io, " %6u\n", length(kpoints))
		for k in kpoints
			@printf(io, " %20.12f %20.12f %20.12f 1\n", k...)
		end
	end
	println(io, "")

	return nothing
end
