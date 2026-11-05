export qₑ, mₑ, mₚ, ħ, ϵ₀, v_light, μB, kB, a₀, Ha, Ry, v_au, t_au, α, CoulombScale

# type declarations on global variables are not yet supported for 1.6.7
"""
	qₑ = 1.602176634

Elementary charge, stored in units of ``10^{-19}`` C.
"""
const qₑ = 1.602176634 #e-19 C
"""
	mₑ = 9.1093837015

Electron rest mass, stored in units of ``10^{-31}`` kg.
"""
const mₑ = 9.1093837015 #e-31 kg
"""
	mₚ = 1.67262192369

Proton rest mass, stored in units of ``10^{-27}`` kg.
"""
const mₚ = 1.67262192369 #e-27 kg
"""
	ħ = 1.054571817646

Reduced Planck constant, stored in units of ``10^{-34}`` J⋅s.
"""
const ħ = 1.054571817646 #e-34 J⋅s
"""
	ϵ₀ = 8.854187817

Vacuum permittivity, stored in units of ``10^{-12}`` F/m.
"""
const ϵ₀ = 8.854187817 #e-12 F/m
"""
	v_light = 299792458

Speed of light in vacuum, in m/s.
"""
const v_light = 299792458 #m/s
"""
	μB = 9.2740100783

Bohr magneton, stored in units of ``10^{-24}`` J/T.
"""
const μB = 9.2740100783 #e-24 J/T
"""
	kB = 1.380649

Boltzmann constant, stored in units of ``10^{-23}`` J/K.
"""
const kB = 1.380649 #e-23 J/K

# atomic unit
"""
	a₀ = 0.5291772111569444

Bohr radius, in Å.
"""
const a₀ = 4π * ϵ₀ * ħ^2 / (10 * mₑ * qₑ^2) # 0.529177Å
"""
	Ha = 4.359744718037628

Hartree energy, stored in units of ``10^{-18}`` J (≈ 27.211 eV).
"""
const Ha = 10 * ħ^2 / (mₑ * a₀^2) # 4.3597 * 1e-18J
const Ha_eV = 10 * Ha / qₑ # 27.211eV
const Ry_eV = Ha_eV / 2 # 13.6eV
"""
	v_au = 2.187691262599066

Atomic unit of velocity, stored in units of ``10^{6}`` m/s.
"""
const v_au = 10 * ħ / (mₑ * a₀) # 2.1877 * 1e6m/s
"""
	t_au = 2.418884328898678

Atomic unit of time, stored in units of ``10^{-17}`` s.
"""
const t_au = mₑ * a₀^2 / ħ # 2.4189 * 1e-17s

"""
	α = 0.007297352533285885

Fine-structure constant (dimensionless).
"""
const α = 1 / 137.03599976 # qₑ²/4πϵ₀cħ

"""
	CoulombScale = e^2/4πϵ₀

Bare Coulomb potential is CoulombScale/r, the unit of r is Å, of potential energy is eV.
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
const _neighbour_offsets_2d = SVector{2,Int}[
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
