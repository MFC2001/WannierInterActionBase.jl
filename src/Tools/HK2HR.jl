"""
	HK2HR
"""
function HK2HR(Hk::AbstractArray{<:Number}, kgrid::RedKgrid; kwargs...)

    nk = size(Hk, 3)
    Hk = [view(Hk, :, :, ik) for ik in 1:nk]

    return HK2HR(Hk, kgrid; kwargs...)
end
function HK2HR(Hk::AbstractVector{<:AbstractMatrix{<:Number}}, kgrid::RedKgrid, lattice::Lattice;
    grid_atol::Real=sqrt(eps(Float64)),
)
    @assert length(Hk) == length(kgrid) "Mismatched Hk and kgrid!"

    R, ndegen = gridindex(Val(:WignerSeitz), lattice, kgrid.kgrid_size; atol=grid_atol)

    nk = length(kgrid)
    norb = size(Hk[1], 1)
    nR = length(R)
    vals = zeros(ComplexF64, norb, norb, nR)
    Threads.@threads for iR in eachindex(R)
        RR = R[iR]
        for ik in 1:nk
            vals[:, :, iR] .+= Hk[ik] .* cispi(-2 * (kgrid[ik] ⋅ RR))
        end
        vals[:, :, iR] ./= nk
    end

    vals = reshape(vals, :)
    is = repeat(1:norb, norb * nR)
    js = repeat(1:norb, inner=norb, outer=nR)
    iRs = repeat(1:nR, inner=norb^2)

    return HR(is, js, iRs, vals, norb, R; ndegen, storage=:full, inbound=true, dofilter=false)
end
# function HK2HR(Hk_irred::AbstractVector{<:AbstractMatrix{<:Number}}, irredkgrid::IrredKgrid, heps::Real = 1e-6)

# 	Nk = length(irredkgrid.redkcoord)
# 	norb = size(Hk_irred[1], 1)
# 	gridindex = map(k -> Vec3{Int}(k .* irredkgrid.kgrid_size), irredkgrid.redkdirect)

# 	Hk = Vector{eltype(Hk_irred)}(undef, Nk)
# 	Threads.@threads for k in eachindex(Hk)
# 		U = irredkgrid.lattsymop[k].U
# 		Hk[k] = Hermitian(U * Hk_irred[irredkgrid.irmap[k]] * U')
# 	end


# 	heps2 = heps^2

# 	path = Matrix{Int}(undef, 5, 0)
# 	value = Vector{ComplexF64}(undef, 0)

# 	lk = ReentrantLock()
# 	Threads.@threads for aimpath in gridindex

# 		t = zeros(ComplexF64, norb, norb)
# 		for k in 1:Nk
# 			t .+= Hk[k][:, :] .* cispi(-2 * (irredkgrid.redkcoord[k] ⋅ aimpath))
# 		end
# 		t ./= Nk

# 		I = abs2.(t) .> heps2
# 		value_t = t[I]

# 		I = Tuple.(CartesianIndices(t)[I])
# 		nt = length(I)
# 		path_t = Matrix{Int}(undef, 5, nt)
# 		for (i, I) in enumerate(I)
# 			path_t[:, i] = [aimpath..., I[1], I[2]]
# 		end

# 		lock(lk) do
# 			path = [path path_t]
# 			append!(value, value_t)
# 		end
# 	end

# 	return HR(path, value; hrsort = true)
# end
