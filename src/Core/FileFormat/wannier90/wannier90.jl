export wannier90_chk, wannier90_eig, wannier90_xsf, wannier90_UNK
export wannier90_amn, wannier90_mmn, wannier90_win
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
	wannier_centres::Vector{Vec3{Float64}}
	wannier_spreads::Vector{Float64}
end
"""
wannier90.eig

	read(path/io, ::Type{wannier90_eig})
	write(path/io, band, ::Type{wannier90_eig}; bandindex::AbstractVector{<:Integer})

- `band` can be a `Eigen`, or a vector of `Eigen` or a VecorMat of real numbers.
- `bandindex` is the index of bands to be written, default is all bands.
"""
struct wannier90_eig <: FileFormat end
"""
wannier90_w.xsf

	read(path/io, ::Type{wannier90_xsf}) -> wannier90_xsf

"""
struct wannier90_xsf <: FileFormat
	cell::Cell{ReducedCoordinates{Float64}}
	origin::Vec3{Float64}
	grid_vecs::Mat3{Float64}
	grid_size::NTuple{3, Int}
	value::Array{Float64, 3}
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
- `nnkpts` is a 5xN matrix, each column is (ik1, ib2, G1, G2, G3) for one k-b pair, it must be provided.
"""
struct wannier90_mmn <: FileFormat end
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

include("./hr.jl")
include("./centres.jl")
