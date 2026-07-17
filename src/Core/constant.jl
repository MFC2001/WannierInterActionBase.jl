export bohr_radius, HartreeEnergy, RydbergEnergy, mₑ, mₚ, qₑ, ϵ₀, ħ, μB, kB, CoulombScale

# type declarations on global variables are not yet supported for 1.6.7
const bohr_radius = 0.529177210903 #Å
const HartreeEnergy = 27.211386245988 #eV
const RydbergEnergy = 13.605693122994 #eV
const mₑ = 9.1093837015 #e-31 kg
const mₚ = 1.67262192369 #e-27 kg
const qₑ = 1.602176634 #e-19 C
const ϵ₀ = 8.854187817 #e-12 F/m
const ħ = 1.054571817646 #e-34 J⋅s
const μB = 9.2740100783 #e-24 J/T
const kB = 1.380649 #e-23 J/K

"""
	CoulombScale = (qₑ * 1000.0) / (4π * ϵ₀)

This term equal to e^2/4πϵ₀ with the unit of r is Å.
Bare Coulomb potential is CoulombScale/r, the unit of r is Å, potential Energy unit is eV.
"""
const CoulombScale = (qₑ * 1000.0) / (4π * ϵ₀)

const _neighbour_offsets_3d = ReducedCoordinates{Int}[
	[-1, -1, -1],
	[-1, -1, 0],
	[-1, -1, 1],
	[-1, 0, -1],
	[-1, 0, 0],
	[-1, 0, 1],
	[-1, 1, -1],
	[-1, 1, 0],
	[-1, 1, 1],
	[0, -1, -1],
	[0, -1, 0],
	[0, -1, 1],
	[0, 0, -1],
	[0, 0, 1],
	[0, 1, -1],
	[0, 1, 0],
	[0, 1, 1],
	[1, -1, -1],
	[1, -1, 0],
	[1, -1, 1],
	[1, 0, -1],
	[1, 0, 0],
	[1, 0, 1],
	[1, 1, -1],
	[1, 1, 0],
	[1, 1, 1],
]
const _neighbour_offsets_2d_xy = ReducedCoordinates{Int}[
	[-1, -1, 0],
	[-1, 0, 0],
	[-1, 1, 0],
	[0, -1, 0],
	[0, 1, 0],
	[1, -1, 0],
	[1, 0, 0],
	[1, 1, 0],
]
const _neighbour_offsets_2d_xz = ReducedCoordinates{Int}[
	[-1, 0, -1],
	[-1, 0, 0],
	[-1, 0, 1],
	[0, 0, -1],
	[0, 0, 1],
	[1, 0, -1],
	[1, 0, 0],
	[1, 0, 1],
]
const _neighbour_offsets_2d_yz = ReducedCoordinates{Int}[
	[0, -1, -1],
	[0, -1, 0],
	[0, -1, 1],
	[0, 0, -1],
	[0, 0, 1],
	[0, 1, -1],
	[0, 1, 0],
	[0, 1, 1],
]
const _neighbour_offsets_2d = SVector{2, Int}[
	[-1, -1],
	[-1, 0],
	[-1, 1],
	[0, -1],
	[0, 1],
	[1, -1],
	[1, 0],
	[1, 1],
]
const _neighbour_offsets_1d_x = ReducedCoordinates{Int}[
	[-1, 0, 0],
	[1, 0, 0],
]
const _neighbour_offsets_1d_y = ReducedCoordinates{Int}[
	[0, -1, 0],
	[0, 1, 0],
]
const _neighbour_offsets_1d_z = ReducedCoordinates{Int}[
	[0, 0, -1],
	[0, 0, 1],
]
const _neighbour_offsets_1d = Int[
	-1,
	1,
]
