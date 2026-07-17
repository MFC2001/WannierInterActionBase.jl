
"""
	write(path::AbstractString, data, ::Type{T}; mode = "w", kwargs...) where {T <: FileFormat}
"""
function Base.write(path::AbstractString, data, ::Type{T}; mode = "w", kwargs...) where {T <: FileFormatIO}
	mkpath(dirname(path))
	return open(path, mode) do io
		write(io, data, T; kwargs...)
	end
end

include("./wannier90/wannier90.jl")
include("./BAND_dat.jl")
include("./dat_eigenvalue.jl")
include("./PBAND_dat.jl")
include("./POSCAR.jl")
include("./QEpw.jl")
include("./WanIntijkl.jl")
