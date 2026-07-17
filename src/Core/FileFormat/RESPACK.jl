"""
The format of output file `dir-intW/dat.Wmat` of RESPACK, also `dir-intW/dat.Vmat`.

	read(path/io, ::Type{RESPACKU}; norb::Integer) -> HR

- `norb::Integer`, the number of orbitals.

Compared to [`RESPACKJ`](@ref), this method will read all the data in the file including onsite term.
"""
struct RESPACKU <: FileFormat end
"""
The format of output file `dir-intJ/dat.Jmat` of RESPACK, also `dir-intW/dat.Xmat`.

	read(path::AbstractString/io::IO, ::Type{RESPACKJ}; norb::Integer) -> HR

- `norb`::Integer, the number of orbitals.

Compared to [`RESPACKU`](@ref), this method will read all the data in the file except onsite term.
"""
struct RESPACKJ <: FileFormat end
"""
File `dir-wfn/dat.eigenvalue` contains electronic band structure, output by wan2respack or RESPACK.

	read(path/io, ::Type{dat_eigenvalue}) -> Matrix{Float64}
	write(path/io, band::Matrix{<:Real}, ::Type{dat_eigenvalue})

"""
struct dat_eigenvalue <: FileFormat end
