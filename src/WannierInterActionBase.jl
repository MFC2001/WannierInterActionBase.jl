module WannierInterActionBase

using LinearAlgebra
using SpecialFunctions
using Dates: Dates
using Requires

#core
using StaticArrays
using StructEquality

#algrism
using QuadGK
using QuadOsc
using HCubature
using Optim
using ADTypes

#io
using DelimitedFiles
using Printf

using HDF5
using JLD2
using NCDatasets

using FortranFiles
using XML: XML

include("./Core/Core.jl")
include("./AbstractInterAction/AbstractInterAction.jl")
include("./Dielectric/Dielectric.jl")
include("./IdealLattice/IdealLattice.jl")
include("./wannier/wannier.jl")
include("./IO/IO.jl")
include("./Tools/Tools.jl")

function __init__()
	@require FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341" include("../ext/FFTW/FFTW.jl")
	@require Struve = "bbbc800e-c7bf-11e8-107b-9b6d1d23829a" include("../ext/Struve.jl")
	@require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("../ext/Plots/Plots.jl")
end

end
