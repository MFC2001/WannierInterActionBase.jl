export WannierPWKgrid, translate!
struct WannierPWKgrid{WT <: Number} # PWKgridStates{WT <: Number}
	nk::Int
	nw::Int
	nspin::Int
	nG::Vector{Int} # length.(Gidx)
	nG_max::Int # maximum(nG)
	kgrid::RedKgrid
	Gidx::Vector{Vector{Int}} #[nk][nG]
	Glist::Vector{ReducedCoordinates{Int}}
	components::Vector{Matrix{WT}} # [nk][nG, nw] or [nk][nG*2, nw]
	fft_grid::Ref{NTuple{3, Int}}
	Gmap::Vector{Int}
	period::Vec3{Bool}
	lattice::Lattice{Float64}
	rlattice::ReciprocalLattice{Float64}
	ldot::Mat3{Float64}
	rldot::Mat3{Float64}
	Ω::Float64
	rΩ::Float64
end
function Base.getproperty(wave::WannierPWKgrid{WT}, name::Symbol) where {WT}
	if name === :fft_grid
		# 解引用 Ref，直接返回 NTuple
		return getfield(wave, :fft_grid)[]
	else
		return getfield(wave, name)
	end
end
function Base.setproperty!(wave::WannierPWKgrid{WT}, name::Symbol, value) where {WT}
	if name === :fft_grid
		value = (value[1], value[2], value[3])
		getfield(wave, :fft_grid)[] = value
		pgi = PeriodicGridIndex(value)
		wave.Gmap .= pgi.(Val(:linear), wave.Glist)
		return value
	else
		error("Only allow change the fft_grid of WannierPWKgrid!")
	end
end
function WannierPWKgrid(kgrid::Union{MonkhorstPack, RedKgrid}, components::Vector{Matrix{WT}},
	Gidx::Vector{Vector{Int}}, Glist::Vector{ReducedCoordinates{Int}}, lattice::Lattice;
	nspin::Integer = 1, fft_grid = (0, 0, 0), period = [true, true, true],
) where {WT <: Number}
	check_lattice_period(lattice, period)
	# 对于 nspin=4 的旋量波函数，前半和后半为不同自旋。
	if kgrid isa MonkhorstPack
		kgrid = RedKgrid(kgrid)
	end
	nk = length(kgrid)
	nw = size(components[1], 2)
	nG = length.(Gidx)
	nG_max = maximum(nG)
	fft_grid = (fft_grid[1], fft_grid[2], fft_grid[3])
	if all(iszero, fft_grid)
		Gmap = Vector{Int}(undef, length(Glist))
	else
		pgi = PeriodicGridIndex(fft_grid)
		Gmap = pgi.(Val(:linear), Glist)
	end
	period = Vec3{Bool}(collect(period))
	# lattice
	lattice = convert(Lattice{Float64}, lattice)
	rlattice = reciprocal(lattice)
	ldot = parent(lattice)
	ldot = transpose(ldot) * ldot
	rldot = parent(rlattice)
	rldot = transpose(rldot) * rldot
	Ω = abs(det(parent(lattice)))
	rΩ = abs(det(parent(rlattice)))
	return WannierPWKgrid{WT}(
		nk, nw, nspin, nG, nG_max, kgrid, Gidx, Glist, components, fft_grid, Gmap,
		period, lattice, rlattice, ldot, rldot, Ω, rΩ,
	)
end
(wave::WannierPWKgrid)(sym::Symbol, args...; kwargs...) = wave(Val(sym), args...; kwargs...)
"""
	(wave::WannierPWKgrid{WT})(::Val{:spin_polarization}) -> Vector{real(WT)}

Return ψ_up^2 - ψ_dn^2, which is the expectation value of σz.
"""
function (wave::WannierPWKgrid{WT})(::Val{:spin_polarization}) where {WT}
	@assert wave.nspin == 4 "Only spinor wave function need to analyse spin structure!"
	ψ_spin_polarization = zeros(real(WT), wave.nw)
	for ik in Base.OneTo(wave.nk)
		nG_ik = wave.nG[ik]
		ψ_ik = wave.components[ik]
		for iw in Base.OneTo(wave.nw)
			ψ_iw_ik_up = view(ψ_ik, 1:nG_ik, iw)
			ψ_iw_ik_dn = view(ψ_ik, (nG_ik+1):(2*nG_ik), iw)
			ψ_spin_polarization[iw] += (sum(abs2, ψ_iw_ik_up) - sum(abs2, ψ_iw_ik_dn))
		end
	end
	ψ_spin_polarization ./= wave.nk
	return ψ_spin_polarization
end
"""
	(wave::WannierPWKgrid)(::Val{:spin_index}; atol::Real = 1e-6) -> (upindex, dnindex)
1 - |ψ_up^2 - ψ_dn^2| < atol, then ψ obtain a spin index.
"""
function (wave::WannierPWKgrid)(::Val{:spin_index}; atol::Real = 1e-6)
	@assert wave.nspin == 4 "Only spinor wave function need to analyse spin structure!"

	ψ_spin_polarization = wave(Val(:spin_polarization))

	upindex = Vector{Int}(undef, 0)
	dnindex = Vector{Int}(undef, 0)
	spinindex = (upindex, dnindex)

	for iw in Base.OneTo(wave.nw)
		if !signbit(ψ_spin_polarization[iw])
			ispin = 1
		else
			ispin = 2
		end
		if (1 - abs(ψ_spin_polarization[iw])) < atol
			push!(spinindex[ispin], iw)
		end
	end
	return spinindex
end
"""
	translate!(wave::WannierPWKgrid, centres::Pair{<:Integer, <:AbstractVector{<:Integer}}...)

You only need to input the wannier centres R which isn't in the [0, 0, 0] cell.
Before, ψ_{i0} ∝ ∑_{k} ψ_{ik}, but ψ_{i-R} is in the [0, 0, 0] cell; 
change the components in `wave` directly, make sure that ∑_{k} ψ′_{ik} is the [0, 0, 0] cell.

See also `translate`.
"""
function translate!(wave::WannierPWKgrid{WT},
	centres::Pair{<:Integer, <:AbstractVector{<:Integer}}...) where {WT}
	for (iw, R) in centres
		for ik in Base.OneTo(wave.nk)
			k = wave.kgrid[ik]
			enkR = cispi(-2 * (k ⋅ R))
			ψ_ik_iw = view(wave.components[ik], :, iw)
			ψ_ik_iw .*= enkR
		end
	end
	return wave
end
function wannier_spin_projected(wave::WannierPWKgrid{WT}; pair::Union{AbstractVector{NTuple{2, <:Integer}}, Nothing} = nothing) where {WT}
	#TODO 理论上，目标子空间存在一组基全为Sz的本征态，充要条件为该子空间为Sz的不变子空间，至少，不是任意体系都可以做到。
	#TODO 这里可以应用optim方法，目标函数可以为每个wannier的(up^2-dn^2)^2，求该函数最大值，可以转化为最小值问题。
	#TODO 为了不破坏局域性，只能提供一个幺正变换矩阵，而不是每个ik一个，也就是等价于对实空间的基进行线性组合。
	#TODO 无需变换至实空间，但是需要明确wannier基属于哪个元胞，可以使用translate函数进行预操作。
	#TODO 参数为每个k点的幺正矩阵，初猜需要用户提供，也就是用户认为可以怎么变换一下，也可以直接选择恒等矩阵。
	return newwave, U
end
