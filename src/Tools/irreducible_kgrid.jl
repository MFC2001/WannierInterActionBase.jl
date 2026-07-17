
function irreducible_kgrid(kgrid::MonkhorstPack, symops::AbstractVector{<:SymOp}, cell::Cell, orbital::wannier90_centres; symtol = 1e-5)

	irredkgrid = irreducible_kgrid(kgrid, symops)

	try
		get_lattsymop!(irredkgrid, orbital.location; symtol)
	catch err
		println(println(err.msg, " Try set a smaller symtol."))
	end

	return irredkgrid
end

function irreducible_kgrid(kgrid::MonkhorstPack, symops::AbstractVector{<:SymOp}, TB::AbstractTightBindModel; symtol = 1e-5)

	irredkgrid = irreducible_kgrid(kgrid, symops)


	get_lattsymop!(irredkgrid, TB.orb_location; symtol)


	return irredkgrid
end

"""
Construct the irreducible wedge given the crystal `symmetries`. Returns the list of k-point
coordinates and the associated weights.
"""
function irreducible_kgrid(kgrid::MonkhorstPack, symops::AbstractVector{<:SymOp})

	if all(isone, kgrid.kgrid_size)
		return IrredKgrid([Vec3(0, 0, 0)], [1], [Vec3(0, 0, 0)], [1],
			one(SymOp), LattSymOp{ComplexF64}[], kgrid.kgrid_size, kgrid.kshift)
	end


	allkdirect = RedKgrid(kgrid).kdirect
	irredk = similar(allkdirect, 0)

	Nk = length(allkdirect)
	irmap = zeros(Int, Nk)
	ksign = ones(Bool, Nk)
	ksymop = Vector{SymOp}(undef, Nk)

	n = 0
	while count(ksign) > 0
		i = findfirst(ksign)
		push!(irredk, allkdirect[i])
		n += 1
		ksign[i] = false
		irmap[i] = n
		ksymop[i] = one(SymOp)

		for symop in symops
			k′ = symop.S * allkdirect[i]
			I = findfirst(allkdirect) do kk
				all(isinteger, kk - k′)
			end
			if isnothing(I) || !ksign[I]
				continue
			else
				ksign[I] = false
				irmap[I] = n
				ksymop[I] = symop
			end
		end
	end

	kweight = [count(irmap .== i) // Nk for i in eachindex(irredk)]

	return IrredKgrid(irredk, kweight, allkdirect, irmap,
		ksymop, similar(ksymop, LattSymOp{ComplexF64}), kgrid.kgrid_size, kgrid.kshift)
end


function get_lattsymop!(irredkgrid::IrredKgrid, orblocation; symtol = 1e-5)
	# The orblocation should be reduced coordinates.
	irredkgrid.lattsymop .= map(zip(irredkgrid.symop, irredkgrid.irmap, irredkgrid.redkdirect)) do (symop, irredk, redk)
		return symop2lattsymop(symop, orblocation, irredkgrid.kdirect[irredk], redk; symtol)
	end
	# irredkgrid.alllattsymop .= map(CartesianIndices(irredkgrid.symmap)) do I
	# 	(i, j) = Tuple(I)
	# 	return symop2lattsymop(irredkgrid.allsymop[i], orblocation, TBmodel.basis, irredkgrid.kcoord[j]; symtol)
	# end
	return irredkgrid
end



function unfold_kcoords(kcoords, symmetries; symtol = 1e-5)
	# unfold
	all_kcoords = [normalize_kpoint_coordinate(symop.S * kcoord)
				   for kcoord in kcoords, symop in symmetries]
	# uniquify
	digits = ceil(Int, -log10(symtol))
	return unique(all_kcoords) do k
		# if x and y are both close to a round value, round(x)===round(y), except at zero
		# where 0.0 and -0.0 are considered different by unique. Add 0.0 to make both
		# -0.0 and 0.0 equal to 0.0
		normalize_kpoint_coordinate(round.(k; digits) .+ 0.0)
	end
end
