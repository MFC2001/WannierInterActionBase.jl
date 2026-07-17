export wannier90_nnkp, wannier90_chk,
	wannier90_eig, wannier90_xsf, wannier90_UNK, wannier90_amn, wannier90_mmn, wannier90_win, wannier90_hr
"""
wannier90.chk

	read(path/io, ::Type{wannier90_chk}) -> wannier90_chk

"""
struct wannier90_chk <: FileFormat
	header::String
	num_bands::Int
	exclude_bands::Vector{Int}
	include_bands::Vector{Int}
	inwindow_bands::Vector{Vector{Int}}
	lattice::Lattice{Float64}
	rlattice::ReciprocalLattice{Float64}
	mp_grid::Vector{Int}
	kdirect::Vector{Vec3{Float64}}
	nntot::Int
	num_wann::Int
	checkpoint::String
	have_disentangled::Bool
	omega_invariant::Float64
	lwindow::Matrix{Bool}
	ndimwin::Vector{Int}
	u_matrix_opt::Vector{Matrix{ComplexF64}}
	u_matrix::Vector{Matrix{ComplexF64}}
	m_matrix::Array{ComplexF64, 4}
	wannier_centres::Vector{ReducedCoordinates{Float64}}
	wannier_spreads::Vector{Float64}
end
"""
wannier90.nnkp

	read(path/io, ::Type{wannier90_nnkp}) -> wannier90_nnkp

"""
struct wannier90_nnkp
	lattice::Lattice{Float64}
	rlattice::ReciprocalLattice{Float64}
	kpoints::Vector{ReducedCoordinates{Float64}}
	nnkpts::Matrix{Int}
	exclude_bands::Vector{Int}
end
"""
wannier90.eig

	read(path/io, ::Type{wannier90_eig}) -> Matrix{Float64}

- return a Matrix with size of (nband, nk).

	write(path/io, band, ::Type{wannier90_eig}; bandindex::AbstractVector{<:Integer})

- `band` can be a `Eigen`, or a vector of `Eigen` or a VecorMat of real numbers.
- `bandindex` is the index of bands to be written, default is all bands.
"""
struct wannier90_eig <: FileFormat end
"""
wannier90_w.xsf

	read(path/io, ::Type{wannier90_xsf}; period = [1, 1, 1]) -> wannier90_xsf{Float64}
	write(path/io, wannier::wannier90_xsf{<:Number}, ::Type{wannier90_xsf}; value = real)

- `period` is the period of unitcell, is used to construct a cell.
- `value` is a function that can be applied to the value of `wannier`, default is `real`.
"""
struct wannier90_xsf{T <: Number} <: FileFormat
	cell::Cell{ReducedCoordinates{Float64}}
	origin::Vec3{Float64}
	grid_vecs::Mat3{Float64}
	grid_size::NTuple{3, Int}
	value::Array{T, 3}
end
"""
	wannier90_xsf(value::AbstractArray{T, 3}, range, cell::Cell{P}) -> wannier90_xsf{T}

- `range` contains three integers, which are the number of unitcells in three directions.
"""
function wannier90_xsf(value::AbstractArray{T, 3}, range, cell::Cell{P}) where {T <: Number, P}
	lattice = parent(cell.lattice)
	origin = -lattice * [floor.(Int, (range .- 1) ./ 2)...]
	origin = Vec3(origin)

	grid_size = size(value)
	grid_vecs = Matrix{Float64}(undef, 3, 3)
	for i in 1:3
		grid_vecs[:, i] = lattice[:, i] * range[i] * (grid_size[i] - 1) / grid_size[i]
	end

	if P <: CartesianCoordinates
		cell = convert(cell)
	end

	return wannier90_xsf{T}(cell, origin, grid_vecs, grid_size, value)
end
"""
UNK_ik.ispin

	read(path/io, ::Type{wannier90_UNK}) -> wannier90_UNK

"""
struct wannier90_UNK <: FileFormat
	ik::Int
	nbnd::Int
	fft_grid::NTuple{3, Int}
	u::Vector{Array{ComplexF64, 3}}
end
"""
wannier90.amn

	write(path/io, A::AbstractArray{<:Number}, ::Type{wannier90_amn})		

(nband, nwannier, nk) = size(A)
"""
struct wannier90_amn <: FileFormat end
"""
wannier90.mmn

	write(path/io, M::AbstractArray{<:Number}, ::Type{wannier90_mmn}; nnkpts::AbstractMatrix{<:Integer})

- `M` is a 3D array, the first two dimensions are band indices, the third dimension is k-b pair index.
- `nnkpts` is a 5xN matrix, each column is (ik, ikb, G1, G2, G3) for one k-b pair, it must be provided.
"""
#TODO
struct wannier90_mmn{MT <: Number} <: FileFormat
	M::Array{MT, 4} # [nb, nb, nntot, nk]
	nnkpts::Matrix{Int} # [5, nntot * nk], [ik, ikb, G1, G2, G3] 同个ik连续存储
	blists::Matrix{CartesianCoordinates{Float64}} # [nntot, nk]
	weights::Matrix{Float64} # [nntot, nk]
end
"""
wannier90.win

	write(io::IO, cell::Cell, ::Type{wannier90_win}; kwargs...)

kwargs:
- `num_iter`: ;
- `num_print_cycles`: default is 20;
- `search_shells`: default is 100;
- `shell_list`: default is `nothing`;
- `dis_num_iter`: ;
- `dis_conv_tol`: default is 1e-12;
- `dis_conv_window`: default is 5;
- `exclude_bands`: default is `nothing`;
- `num_bands`: ;
- `dis_win_min`: ;
- `dis_win_max`: ;
- `dis_froz_min`: ;
- `dis_froz_max`: ;
- `num_wann`: ;
- `use_bloch_phases`: default is `false`;
- `spinors`: default is `false`;
- `guiding_centres`: default is `false`;
- `projections`: default is `nothing`;
- `wannier_plot`: default is `false`;
- `wannier_plot_format`: default is "xcrysden";
- `wannier_plot_supercell`: default is `nothing`;
- `bands_plot`: default is `false`;
- `bands_num_points`: default is 30;
- `bands_plot_format`: default is "gnuplot xmgrace";
- `kpoint_path`: default is `nothing`, its format is the same as input of [`Kline`](@ref);
- `grid`: Union{[`MonkhorstPack`](@ref), [`RedKgrid`](@ref)}.
"""
struct wannier90_win <: FileFormat end
"""
wannier90_hr.dat

	read(path/io, ::Type{wannier90_hr}; atol::Real = 0, readimag::Bool = true, μ::Real = 0, hrsort::Bool = false) -> HR

- `atol`: can be set as a random positive real number, which is used to filter small values.
- `readimag`: decides whether to read imaginary part.
- `μ`: can be set as a random real number, which is the chemical potential.
- `hrsort`: decides whether to sort the object according to the standard order.

	write(path/io, hr::HR, ::Type{wannier90_hr})

Write hr to `path/io` in the format of `wannier90_hr.dat`.
"""
struct wannier90_hr <: FileFormat end


include("./centres.jl")
