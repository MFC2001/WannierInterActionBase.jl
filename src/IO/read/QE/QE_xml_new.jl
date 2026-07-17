function Base.read(path::AbstractString, ::Type{QE_xml_new}; period = [1, 1, 1])

	xml = read(path, XML.Node)
	xml = xml[end]

	data = Dict{String, Any}()

	atomic_structure = xml["output"][1]["atomic_structure"][1]
	nat = parse(Int, XML.attributes(atomic_structure)["nat"])
	alat = parse(Float64, XML.attributes(atomic_structure)["alat"]) * bohr_radius
	ralat = 2π / alat
	lattice = atomic_structure["cell"][1]
	a1 = parse.(Float64, split(XML.value(XML.children(lattice["a1"][1])[1])))
	a2 = parse.(Float64, split(XML.value(XML.children(lattice["a2"][1])[1])))
	a3 = parse.(Float64, split(XML.value(XML.children(lattice["a3"][1])[1])))
	lattice_data = [a1 a2 a3]
	lattice = Lattice(lattice_data * bohr_radius)
	rlattice = reciprocal(lattice)
	atomic_positions = atomic_structure["atomic_positions"][1]
	location = Vector{ReducedCoordinates{Float64}}(undef, nat)
	name = Vector{String}(undef, nat)
	index = Vector{Int}(undef, nat)
	for (i, atom) in enumerate(atomic_positions["atom"])
		name[i] = XML.attributes(atom)["name"]
		index[i] = parse(Int, XML.attributes(atom)["index"])
		location[i] = ReducedCoordinates(lattice_data \ parse.(Float64, split(XML.value(XML.children(atom)[1]))))
	end
	data["cell"] = Cell(lattice, location, "ReducedCoordinates"; name, index, period)

	# symmetries
	symmetries = xml["output"][1]["symmetries"][1]
	data["nsym"] = parse(Int, XML.value(XML.children(symmetries["nsym"][1])[1]))
	data["nrot"] = parse(Int, XML.value(XML.children(symmetries["nrot"][1])[1]))
	data["space_group"] = XML.value(XML.children(symmetries["space_group"][1])[1])
	symmetries = symmetries["symmetry"]
	nsym = length(symmetries)
	symmetries_info = Vector{String}(undef, nsym)
	symmetries_name = Vector{String}(undef, nsym)
	symmetries_class = Vector{String}(undef, nsym)
	symmetries_rank = Vector{Int}(undef, nsym)
	symmetries_rotation = Vector{Mat3{Int}}(undef, nsym)
	for (i, sym) in enumerate(symmetries)
		info = sym["info"][1]
		symmetries_info[i] = XML.value(XML.children(info)[1])
		info_attributes = XML.attributes(info)
		symmetries_name[i] = info_attributes["name"]
		symmetries_class[i] = get(info_attributes, "class", "")
		rotation = sym["rotation"][1]
		symmetries_rank[i] = parse(Int, XML.attributes(rotation)["rank"])
		rotation = XML.value(XML.children(rotation)[1])
		# 这里的转置可以通过作用在原子分数坐标下来验证。
		symmetries_rotation[i] = reshape(parse.(Float64, split(rotation)), 3, 3)'
	end
	crystal_symmetry = findall(==("crystal_symmetry"), symmetries_info)
	lattice_symmetry = findall(==("lattice_symmetry"), symmetries_info)
	data["symmetries_name"] = symmetries_name[crystal_symmetry]
	data["symmetries_class"] = symmetries_class[crystal_symmetry]
	data["symmetries_rank"] = symmetries_rank[crystal_symmetry]
	crystal_symmetries_rotation = symmetries_rotation[crystal_symmetry]
	data["lattice_symmetries_name"] = symmetries_name[lattice_symmetry]
	# data["lattice_symmetries_class"] = symmetries_class[lattice_symmetry]
	data["lattice_symmetries_rank"] = symmetries_rank[lattice_symmetry]
	data["lattice_symmetries_rotation"] = symmetries_rotation[lattice_symmetry]
	symmetries_mat = Vector{SymOp}(undef, length(crystal_symmetry))
	for (i, isym) in enumerate(crystal_symmetry)
		fractran = XML.value(XML.children(symmetries[isym]["fractional_translation"][1])[1])
		symmetries_mat[i] = SymOp(crystal_symmetries_rotation[i], parse.(Float64, split(fractran)))
	end
	data["symmetries"] = symmetries_mat

	fft_grid = XML.attributes(xml["output"][1]["basis_set"][1]["fft_grid"][1])
	data["fft_grid"] = parse.(Int, (fft_grid["nr1"], fft_grid["nr2"], fft_grid["nr3"]))

	band_structure = xml["output"][1]["band_structure"][1]
	# obtain nspin from lsda, noncolin and spinorbit.
	lsda = parse(Bool, XML.value(XML.children(band_structure["lsda"][1])[1]))
	noncolin = parse(Bool, XML.value(XML.children(band_structure["noncolin"][1])[1]))
	if lsda
		data["nspin"] = 2
	else
		if noncolin
			data["nspin"] = 4
		else
			data["nspin"] = 1
		end
	end
	data["spinorbit"] = parse(Bool, XML.value(XML.children(band_structure["spinorbit"][1])[1]))

	if data["nspin"] == 2
		data["nbnd_up"] = parse(Int, XML.value(XML.children(band_structure["nbnd_up"][1])[1]))
		data["nbnd_dw"] = parse(Int, XML.value(XML.children(band_structure["nbnd_dw"][1])[1]))
		data["nbnd"] = data["nbnd_up"] + data["nbnd_dw"]
	else
		data["nbnd"] = parse(Int, XML.value(XML.children(band_structure["nbnd"][1])[1]))
	end

	data["nelec"] = Int(parse(Float64, XML.value(XML.children(band_structure["nelec"][1])[1])))
	data["occupations_kind"] = XML.value(XML.children(band_structure["occupations_kind"][1])[1])

	if data["nspin"] == 2 && data["occupations_kind"] == "fixed"
		data["two_fermi_energies"] = parse.(Float64, split(XML.value(XML.children(band_structure["two_fermi_energies"][1])[1])))
	else
		data["fermi_energy"] = parse(Float64, XML.value(XML.children(band_structure["fermi_energy"][1])[1]))
	end

	data["nk"] = parse(Int, XML.value(XML.children(band_structure["nks"][1])[1]))


	ks_energies = band_structure["ks_energies"]
	@assert length(ks_energies) == data["nk"]
	kpoints = Vector{ReducedCoordinates{Float64}}(undef, data["nk"])
	kweights = Vector{Float64}(undef, data["nk"])
	npws = Vector{Int}(undef, data["nk"])
	eigenvalues = Matrix{Float64}(undef, data["nbnd"], data["nk"])
	occupations = Matrix{Float64}(undef, data["nbnd"], data["nk"])
	for (ik, ks_energy) in enumerate(ks_energies)
		kpoint = ks_energy["k_point"][1]
		kweights[ik] = parse(Float64, XML.attributes(kpoint)["weight"])
		kpoint = CartesianCoordinates(parse.(Float64, split(XML.value(XML.children(kpoint)[1]))) .* ralat)
		kpoints[ik] = rlattice \ kpoint
		npws[ik] = parse.(Int, XML.value(XML.children(ks_energy["npw"][1])[1]))
		eigenvalues[:, ik] .= parse.(Float64, split(XML.value(XML.children(ks_energy["eigenvalues"][1])[1]))) .* HartreeEnergy
		occupations[:, ik] .= parse.(Float64, split(XML.value(XML.children(ks_energy["occupations"][1])[1])))
	end
	data["kpoints"] = kpoints
	data["kweights"] = kweights
	data["npws"] = npws
	data["eigenvalues"] = eigenvalues
	data["occupations"] = occupations

	if data["occupations_kind"] == "fixed"
		data["occupations"] = Int.(data["occupations"])
	end

	if data["nspin"] == 2
		data["eigenvalues_up"] = data["eigenvalues"][1:data["nbnd_up"], :]
		data["eigenvalues_dn"] = data["eigenvalues"][(data["nbnd_up"]+1):end, :]
	end

	return data
end

@inline function Base.getindex(node::XML.Node, name::String)
	return filter(n -> XML.tag(n) == name, XML.children(node))
end
