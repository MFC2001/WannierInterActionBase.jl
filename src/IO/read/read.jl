
"""
	read(path::AbstractString, ::Type{T}; kwargs...) where {T <: FileFormat}
"""
function Base.read(path::AbstractString, ::Type{T}; kwargs...) where {T <: FileFormat}
	return open(path, "r") do io
		read(io, T; kwargs...)
	end
end

include("./BAND.jl")
include("./POSCAR.jl")
include("./QE_wfc.jl")
include("./QE_xml_new.jl")
include("./RESPACKUJ.jl")
include("./dat_eigenvalue.jl")
include("./wannier90/wannier90.jl")
