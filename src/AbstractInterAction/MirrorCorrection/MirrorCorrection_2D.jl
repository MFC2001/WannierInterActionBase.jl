
function MirrorCorrection_2D(U::HR, kgrid::RedKgrid, lattice::Lattice, orblocat::AbstractVector{<:ReducedCoordinates}, rcut::Real;
	αrcut = 4.5, δ = 1e-8, ϵ = 1, head = nothing, period = Bool[1, 1, 0])

	if length(period) ≠ 3 || count(period) ≠ 2
		error("Wrong period.")
	end

	α = αrcut / rcut

	# Create Ggrid.
	G2_max = -4 * α^2 * log(δ)

	rlattice = reciprocal(lattice)
	b₁ = rlattice[:, 1]
	b₂ = rlattice[:, 2]
	b₃ = rlattice[:, 3]

	V_BZ = abs((b₁ × b₂) ⋅ b₃)

	h₁ = V_BZ / norm(b₂ × b₃)
	h₂ = V_BZ / norm(b₃ × b₁)
	h₃ = V_BZ / norm(b₁ × b₂)
	Ggrid = Int.(cld.(√G2_max * 1.2, [h₁, h₂, h₃])) * 2 .+ 1

	a₁ = lattice[:, 1]
	a₂ = lattice[:, 2]
	a₃ = lattice[:, 3]

	if !period[1]
		Ggrid[1] = 1
		S = norm(a₂ × a₃)
		zindex = 1
	elseif !period[2]
		Ggrid[2] = 1
		S = norm(a₃ × a₁)
		zindex = 2
	elseif !period[3]
		Ggrid[3] = 1
		S = norm(a₁ × a₂)
		zindex = 3
	end

	rmat = parent(rlattice)
	rldot = transpose(rmat) * rmat
	Ggrid = gridindex(Ggrid)
	Ggrid = filter(G -> G ⋅ (rldot * G) < G2_max, Ggrid)

	allkG = reshape([k + G for k in kgrid.kdirect, G in Ggrid], :)

	# gauss potential
	φR = RealGauss(; ϵ, α)
	φK = ReciprocalGauss2D(; ϵ, α, S)

	# pre-calculation
	nk = length(kgrid)
	norb = length(orblocat)
	φkG = Array{Float64}(undef, length(allkG), norb, norb)
	φ₀ = Matrix{Float64}(undef, norb, norb)

	if isnothing(head)
		allkG_car = map(kG -> rlattice * kG, allkG)
		for i in 1:norb, j in 1:i
			dorb = lattice * (orblocat[j] - orblocat[i])
			z = dorb[zindex]
			φkG[:, i, j] = map(kG -> φK(kG, z), allkG_car)
			φkG[:, j, i] = φkG[:, i, j]
			φ₀[i, j] = φK(Val(:head), nk, z) / nk
			φ₀[j, i] = φ₀[i, j]
		end
	else
		φ₀ .= head
	end

	φ_LR(r, i, j) = sum(kGi -> φkG[kGi, i, j] * cospi(2 * (allkG[kGi] ⋅ r)), eachindex(allkG)) / nk + φ₀[i, j]

	correction = Vector{Float64}(undef, length(U.value))
	Threads.@threads for i in eachindex(U.value)
		path = ReducedCoordinates(U.path[1, i], U.path[2, i], U.path[3, i])
		i_idx = U.path[4, i]
		j_idx = U.path[5, i]
		r_frac = path + orblocat[j_idx] - orblocat[i_idx]
		correction[i] = φ_LR(r_frac, i_idx, j_idx) - φR(lattice * r_frac)
	end

	U_corrected = HR(U.path, U.value - correction; hrsort = false)

	return U_corrected
end
