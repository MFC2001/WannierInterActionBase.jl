
function Base.read(io::IO, ::Type{QE_xml_new}; period = [1, 1, 1])
	filename = io.name[2:end-1]
	close(io)
	if filename[1:4] == "file"
		filename = filename[6:end]
	else
		error("Only support XML file.")
	end
	xml = read(filename, XML.Node)
	xml = xml[end]

	atomic_structure = xml["output"][1]["atomic_structure"][1]
	nat = parse(Int, XML.attributes(atomic_structure)["nat"])
	alat = parse(Float64, XML.attributes(atomic_structure)["alat"]) * bohr_radius
	lattice = atomic_structure["cell"][1]
	a1 = parse.(Float64, split(XML.value(XML.children(lattice["a1"][1])[1])))
	a2 = parse.(Float64, split(XML.value(XML.children(lattice["a2"][1])[1])))
	a3 = parse.(Float64, split(XML.value(XML.children(lattice["a3"][1])[1])))
	lattice_data = [a1 a2 a3]
	lattice = Lattice(lattice_data * bohr_radius)
	atomic_positions = atomic_structure["atomic_positions"][1]
	location = Vector{ReducedCoordinates{Float64}}(undef, nat)
	name = Vector{String}(undef, nat)
	index = Vector{Int}(undef, nat)
	for (i, atom) in enumerate(atomic_positions["atom"])
		name[i] = XML.attributes(atom)["name"]
		index[i] = parse(Int, XML.attributes(atom)["index"])
		location[i] = ReducedCoordinates(lattice_data \ parse.(Float64, split(XML.value(XML.children(atom)[1]))))
	end
	cell = Cell(lattice, location, "ReducedCoordinates"; name, index, period)

	fft_grid = XML.attributes(xml["output"][1]["basis_set"][1]["fft_grid"][1])
	fft_grid = parse.(Int, (fft_grid["nr1"], fft_grid["nr2"], fft_grid["nr3"]))

	band_structure = xml["output"][1]["band_structure"][1]
	nk = parse(Int, XML.value(XML.children(band_structure["nks"][1])[1]))
	# ks_energies = band_structure["ks_energies"]
	# kpoints = map(ks_energies) do ks_energie
	# 	kpoint = parse.(Float64, split(XML.value(XML.children(ks_energie["k_point"][1])[1])))
	# 	ReducedCoordinates(kpoint)
	# end

	# 根据lsda，noncolin，spinorbit判断波函数类型
	lsda = parse(Bool, XML.value(XML.children(band_structure["lsda"][1])[1]))
	noncolin = parse(Bool, XML.value(XML.children(band_structure["noncolin"][1])[1]))
	spinorbit = parse(Bool, XML.value(XML.children(band_structure["spinorbit"][1])[1]))
	nspin = nothing
	if lsda
		nspin = 2
	else
		if noncolin
			nspin = 4
		else
			nspin = 1
		end
	end

	return QE_xml_new(cell, fft_grid, nk, nspin)
end

@inline function Base.getindex(node::XML.Node, name::String)
	return filter(n -> XML.tag(n) == name, XML.children(node))
end
