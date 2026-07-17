export γgraphyneCell, γgraphyneHR
function γgraphyneCell(; a_intra = 1, a_inter = 2, c = 20)

	l = 2 * a_intra + a_inter
	𝐚 = Vec3(l, 0, 0)
	𝐛 = l * Vec3(-1 / 2, √3 / 2, 0)
	𝐜 = Vec3(0, 0, c)

	lattice = Lattice([𝐚 𝐛 𝐜])

	location5 = [a_inter / 4, a_inter * √3 / 4, c / 2]
	location6 = location5 + [-a_intra / 2, a_intra * √3 / 2, 0]
	location4 = location5 + [a_intra, 0, 0]

	location3 = location4 + [a_intra / 2, a_intra * √3 / 2, c / 2]
	location2 = location3 + [-a_intra / 2, a_intra * √3 / 2, 0]
	location1 = location6 + [a_intra / 2, a_intra * √3 / 2, 0]

	location = CartesianCoordinates.([location1, location2, location3, location4, location5, location6])
	location = map(x -> lattice \ x, location)

	return Cell(lattice, location; location_type = "ReducedCoordinates", name = fill("C", 6), period = Bool[1, 1, 0])
end
function γgraphyneHR(; t₀ = 0, t_intra = -2.27, t_inter = 1.8)
	allpath = Matrix{Int}(undef, 5, 0)
	allvalue = Vector{Float64}(undef, 0)
	if !iszero(t₀)
		path_t₀ = [
			0 0 0 1 1;
			0 0 0 2 2;
			0 0 0 3 3;
			0 0 0 4 4;
			0 0 0 5 5;
			0 0 0 6 6
		]
		value_t₀ = fill(t₀, size(path_t₀, 1))
		allpath = [allpath transpose(path_t₀)]
		allvalue = [allvalue; value_t₀]
	end
	if !iszero(t_intra)
		path_t_intra = [
			0 0 0 1 2;
			0 0 0 1 6;
			0 0 0 2 1;
			0 0 0 2 3;
			0 0 0 3 2;
			0 0 0 3 4;
			0 0 0 4 3;
			0 0 0 4 5;
			0 0 0 5 4;
			0 0 0 5 6;
			0 0 0 6 5;
			0 0 0 6 1
		]
		value_t_intra = fill(t_intra, size(path_t_intra, 1))
		allpath = [allpath transpose(path_t_intra)]
		allvalue = [allvalue; value_t_intra]
	end
	if !iszero(t_inter)
		path_t_inter = [
			0 1 0 1 4;
			1 1 0 2 5;
			1 0 0 3 6;
			0 -1 0 4 1;
			-1 -1 0 5 2;
			-1 0 0 6 3
		]
		value_t_inter = fill(t_inter, size(path_t_inter, 1))
		allpath = [allpath transpose(path_t_inter)]
		allvalue = [allvalue; value_t_inter]
	end

	return HR(allpath, allvalue; orbindex = collect(1:6), hrsort = true)
end
