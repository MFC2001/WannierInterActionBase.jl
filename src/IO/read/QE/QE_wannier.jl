function Base.read(savefolder::AbstractString, chkdata::wannier90_chk, ::Type{QE_wannier}; kwargs...)
	return read(savefolder, QE_wannier; kwargs...,
		bandindex = chkdata.inwindow_bands, u_matrix = chkdata.u_matrix, u_matrix_opt = chkdata.u_matrix_opt)
end
function Base.read(savefolder::AbstractString, ::Type{QE_wannier};
	wfcspin = nothing, period = [1, 1, 1],
	kgrid::Union{MonkhorstPack, Nothing} = nothing, kgrid_atol::Real = 1e-6,
	bandindex::AbstractVector{<:AbstractVector{Int}},
	u_matrix::AbstractVector{<:AbstractMatrix{<:Number}},
	u_matrix_opt::AbstractVector{<:AbstractMatrix{<:Number}} = Vector{Matrix{ComplexF64}}(undef, 0),
)

	isdir(savefolder) || error("$savefolder isn't found!")

	if isfile(joinpath(savefolder, "data-file-schema.xml"))
		# New style QE output XML file (data-file-schema.xml).
		xmldata = read(joinpath(savefolder, "data-file-schema.xml"), QE_xml_new; period)
	elseif isfile(joinpath(savefolder, "data-file.xml"))
		# Old style QE output XML file (data-file.xml).
		error("Only support qe with new version at present.")
	else
		error("No XML file founded!")
	end

	nspin = xmldata["nspin"]

	wfcname = nothing
	if isnothing(wfcspin)
		wfcname = "wfc"
	elseif wfcspin == :up
		wfcname = "wfcup"
	elseif wfcspin ∈ [:down, :dw, :dn]
		wfcname = "wfcdw"
	else
		error("Wrong wfcspin.")
	end

	if isfile(joinpath(savefolder, wfcname * "1.dat"))
		wfc_suffix = "dat"
		wfc_format = QEwfc_dat
	elseif isfile(joinpath(savefolder, wfcname * "1.hdf5"))
		wfc_suffix = "hdf5"
		wfc_format = QEwfc_hdf5
	else
		error("No wavefunction file founded!")
	end

	nk = xmldata["nk"]
	kpoints = Vector{ReducedCoordinates{Float64}}(undef, nk)
	components = Vector{Matrix{ComplexF64}}(undef, nk) #psi_{ik}
	Glists = Vector{Vector{ReducedCoordinates{Int}}}(undef, nk)
	if isempty(u_matrix_opt)
		for ik in 1:nk
			wfc_ik = read(joinpath(savefolder, wfcname * "$(ik)." * wfc_suffix), wfc_format)
			kpoints[ik] = wfc_ik.xk # the kpoints from .xml is not right.
			components[ik] = wfc_ik.components[:, bandindex[ik]] * u_matrix[ik]
			Glists[ik] = wfc_ik.Glist
		end
	else
		for ik in 1:nk
			wfc_ik = read(joinpath(savefolder, wfcname * "$(ik)." * wfc_suffix), wfc_format)
			kpoints[ik] = wfc_ik.xk
			components[ik] = (wfc_ik.components[:, bandindex[ik]] * u_matrix_opt[ik]) * u_matrix[ik]
			Glists[ik] = wfc_ik.Glist
		end
	end

	# check norm
	for ik in 1:nk, iw in axes(components[ik], 2)
		wave_norm = norm(components[ik][:, iw])
		@assert isapprox(wave_norm, 1; atol = 1e-6) "The norm of state with ik=$ik and iw=$iw isn't 1."
	end

	# kgrid
	if kgrid isa MonkhorstPack
		kpoints_rational = kpoints_Float_to_Rational(kpoints, kgrid)
		kgrid = RedKgrid(kpoints_rational, kgrid.kgrid_size, kgrid.kshift)
	else
		kpoints_rational = kpoints_Float_to_Rational(kpoints; tol = kgrid_atol)
		kgrid_xyz = map(1:3) do i
			unique(map(k -> mod(k[i], 1), kpoints_rational))
		end
		kgrid_size = length.(kgrid_xyz)
		if prod(kgrid_size) ≠ nk
			println("nk: ", nk)
			println("kgrid_size: ", kgrid_size)
			println("kgrid_xyz: ", kgrid_xyz)
			error("kgrid is not a MP grid.")
		end
		kshift = map(1:3) do i
			any(iszero, kgrid_xyz[i]) ? 0 : 1 // 2
		end
		kgrid_mp = RedKgrid(MonkhorstPack(kgrid_size; kshift))
		if any(kpoints_rational) do k
			isnothing(findfirst(k_mp -> all(isinteger, k_mp - k), kgrid_mp))
		end
			error("kpoints are not consistent with Monkhorst-Pack grid.")
		end
		kgrid = RedKgrid(kpoints_rational, kgrid_size, kshift)
	end

	# pw.x的输出，可能.xml文件与wfc.dat文件一一对应；
	# 但open_grid.x的输出一定不是一一对应的！

	# rldot
	lattice = xmldata["cell"].lattice
	rlattice = reciprocal(lattice)
	rldot = parent(rlattice)
	rldot = transpose(rldot) * rldot

	# wavek_normalize!(Glists, kgrid)
	(Glist, Gidx) = merge_Glist(Glists...; rldot, Gidx = true)

	return WannierPWKgrid(kgrid, components, Gidx, Glist, xmldata["cell"].lattice;
		nspin, fft_grid = xmldata["fft_grid"], period = xmldata["cell"].period)
end
