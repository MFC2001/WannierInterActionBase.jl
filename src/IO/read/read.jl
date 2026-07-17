
"""
	read(path::AbstractString, ::Type{T}; kwargs...) where {T <: FileFormat}
"""
function Base.read(path::AbstractString, ::Type{T}; kwargs...) where {T <: FileFormatIO}
	return open(path, "r") do io
		read(io, T; kwargs...)
	end
end

include("./BerkeleyGW/BerkeleyGW.jl")
include("./QE/QE.jl")
include("./wannier90/wannier90.jl")
include("./yambo/yambo.jl")
include("./BAND_dat.jl")
include("./dat_eigenvalue.jl")
include("./POSCAR.jl")
include("./RESPACKUJ.jl")
include("./WanIntijkl.jl")
