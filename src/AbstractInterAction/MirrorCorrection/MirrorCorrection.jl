export MirrorCorrection

include("./_getU_out_of_rcut.jl")
include("./Gauss.jl")
include("./KeldyshGauss.jl")

#TODO:输入head的方法有问题，截断库伦势和高斯库伦势中除发散项外均可能含有常数项，需要额外计入。
#TODO: 或许需要带着相位因子展开。

"""
	MirrorCorrection(sym::Symbol, U, kgrid, lattice, orblocat; outfolder::String = "./", logfile::String = joinpath(outfolder, "MirrorCorrection.log"), kwargs...)

Try to apply a mirror correction to `U`.

- `sym` is the symbol of the mirror correction method, currently supports `:Gauss_optim`, `:Gauss`, `:KeldyshGauss_optim`, `:KeldyshGauss`.
- `U::Union{HR, AbstractReciprocalHoppings}`
- `kgrid::Union{MonkhorstPack, RedKgrid, AbstractVector{<:Integer}, NTuple{3, <:Integer}}`
- `lattice::Union{Lattice, AbstractMatrix{<:Real}}`
"""
function MirrorCorrection(sym::Symbol, U::Union{HR, AbstractReciprocalHoppings},
	kgrid::Union{MonkhorstPack, RedKgrid, AbstractVector{<:Integer}, NTuple{3, <:Integer}},
	lattice::Union{Lattice, AbstractMatrix{<:Real}}, orblocat::AbstractVector{V};
	outfile::Union{Nothing, AbstractString} = nothing,
	logfile::AbstractString = isnothing(outfile) ? "./MirrorCorrection.log" : joinpath(dirname(outfile), "MirrorCorrection.log"),
	kwargs...,
) where {V}

	logio = get_logio(logfile)
	original_stdout = stdout
	redirect_stdout(logio)
	_exit_flush(logio)

	# preprocessing
	if U isa HR
		U = ReciprocalHoppings(U)
	end
	if kgrid isa AbstractVector || kgrid isa NTuple
		@assert all(>(0), kgrid) "kgrid should be positive"
		kgrid = RedKgrid(MonkhorstPack(kgrid[1], kgrid[2], kgrid[3]))
	else
		kgrid = RedKgrid(MonkhorstPack(kgrid.kgrid_size))
	end
	lattice = Lattice(parent(lattice))
	orblocat = collect(orblocat)
	if V <: CartesianCoordinates
		orblocat = map(orblocat) do r
			lattice \ r
		end
	elseif V <: AbstractVector
		orblocat = map(orblocat) do r
			ReducedCoordinates(r)
		end
	else
		error("orblocat should be a vector of vector!")
	end

	# calculate
	result = MirrorCorrection(Val(sym), U, kgrid, lattice, orblocat; kwargs...)

	redirect_stdout(original_stdout)
	close(logio)

	if !isnothing(outfile)
		mkpath(dirname(outfile))
		jldsave(outfile; result...)
	end

	return result
end

# 虽然增加高斯展宽后的库伦势在实空间会很快趋于点电荷的结果，但是转换至有限大周期体系中的形式后，高斯展宽参数会体现出较为明显的影响，这为依据远距离的样本点进行参数优化提供了理论依据。
# 2D的要计算远一点。wannier基在非周期方向存在位移的情况下，KeldyshGauss也可以描述，但是需要更远的距离。
# 对比MoS2的结果，wannier基局域性更好后，Mirrorcorrection更难收敛。
