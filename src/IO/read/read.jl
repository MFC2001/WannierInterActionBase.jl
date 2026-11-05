function Base.read(path::AbstractString, ::Type{T}; kwargs...) where {T <: FileFormatIO}
	return open(path, "r") do io
		read(io, T; kwargs...)
	end
end

include("./BerkeleyGW/BerkeleyGW.jl")
include("./QE/QE.jl")
include("./RESPACK/RESPACK.jl")
include("./wannier90/wannier90.jl")
include("./yambo/yambo.jl")
include("./BAND_dat.jl")
include("./POSCAR.jl")
include("./wannierberri/wannierberri.jl")
