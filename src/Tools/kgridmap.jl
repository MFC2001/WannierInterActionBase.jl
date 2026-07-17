"""
	kgridmap(aimkdirects, kdirects, mode = +) -> Matrix{Int}

Try to get a map matrix between aimkdirects and kdirects.

```julia
julia> M = kgridmap(aimkdirects, kdirects, +);

julia> aimkdirects[M[i, j]] == mode(kdirects[i], kdirects[j])
true

```

We also provide some shortcuts:

	kgridmap(kgrid::MonkhorstPack, mode = +)
	kgridmap(redkgrid::RedKgrid, mode = +)

"""
function kgridmap(kgrid::MonkhorstPack, mode::Function = +)
	iszero(kgrid.kshift) && error("Only work for Γ-centered kgrid!")
	redkgrid = RedKgrid(kgrid)
	return kgridmap(redkgrid, redkgrid, mode)
end
function kgridmap(redkgrid::RedKgrid, mode::Function = +)
	iszero(redkgrid.kshift) && error("Only work for Γ-centered kgrid!")
	return kgridmap(redkgrid, redkgrid, mode)
end
function kgridmap(irredkgrid::IrredKgrid, mode::Function = +)
	error("To be continued.")
	return kgridmap(irredkgrid, irredkgrid, mode)
end
kgridmap(aimkdirects, kdirects, mode::Function = +) = kgridmap(aimkdirects, kdirects, kdirects, mode)
function kgridmap(aimkdirects, kdirects1, kdirects2, mode::Function = +)
	# aimk_to_idx = Dict{ReducedCoordinates{Int}, Int}(k => i for (i, k) in enumerate(aimkdirects))

	nk1 = length(kdirects1)
	nk2 = length(kdirects2)

	mapmat = Matrix{Int}(undef, nk1, nk2)
	Threads.@threads for I in CartesianIndices(mapmat)
		@inbounds begin
			(i, j) = Tuple(I)
			k1 = kdirects1[i]
			k2 = kdirects2[j]
			Tk = mode(k1, k2)

			idx = findfirst(aimkdirects) do k
				all(isinteger, k - Tk)
			end
			if isnothing(idx)
				mapmat[I] = 0
				@warn "cannot find mode(k1, k2) in aimkdirects"
			else
				mapmat[I] = idx
			end
		end
	end
	return mapmat
end
