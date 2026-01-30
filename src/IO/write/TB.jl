"""
	write(pathdir::AbstractString, TB::AbstractTightBindModel; mode = "w", comment = "", format = "Cartesian")

Write the `hr.dat`, `POSCAR`, and `centres.xyz` files for a given `AbstractTightBindModel` in the specified directory.
"""
function Base.write(pathdir::AbstractString, TB::AbstractTightBindModel; mode = "w", comment = "", format = "Cartesian")
	write(joinpath(pathdir, "hr.dat"), HR(TB), wannier90_hr; mode, comment)
	write(joinpath(pathdir, "POSCAR"), Cell(TB), POSCAR; mode, comment, format)
	write(joinpath(pathdir, "centres.xyz"), ORBITAL(TB), wannier90_centres; mode, comment)
end
