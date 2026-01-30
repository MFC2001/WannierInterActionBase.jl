
function Base.read(io::IO, ::Type{wannier90_chk}; convert = "native")

	io = FortranFile(io; convert)

	# 1. 读取头部信息 (33字符)
	header = trimstring(read(io, FString{33}))

	# 2. 读取基本参数
	num_bands = Int(read(io, Int32))
	num_exclude_bands = Int(read(io, Int32))

	# 3. 读取排除能带列表
	exclude_bands = nothing
	if iszero(num_exclude_bands)
		read(io)
		read(io) # is it a bug?
		exclude_bands = Int[]
	else
		exclude_bands = read(io, (Int32, num_exclude_bands))
	end
	ib = 0
	include_bands = Vector{Int}(undef, 0)
	sizehint!(include_bands, num_bands)
	while true
		ib += 1
		ib ∈ exclude_bands && continue
		push!(include_bands, ib)
		length(include_bands) == num_bands && break
	end

	# 4. 读取晶格信息
	lattice = Lattice(read(io, (Float64, 3, 3))')
	rlattice = ReciprocalLattice(read(io, (Float64, 3, 3))')

	# 5. 读取k点信息
	num_kpts = Int(read(io, Int32))
	mp_grid = Int.(read(io, (Int32, 3)))
	kdirect = read(io, (Float64, 3, num_kpts))
	kdirect = map(Vec3{Float64}, eachcol(kdirect))

	# 6. 读取其他参数
	nntot = Int(read(io, Int32))
	num_wann = Int(read(io, Int32))
	checkpoint = trimstring(read(io, FString{20}))
	have_disentangled = (read(io, Int32) ≠ 0)

	# 7. 处理disentangled相关数据
	if have_disentangled
		omega_invariant = read(io, Float64)
		lwindow = (read(io, (Int32, num_bands, num_kpts)) .≠ 0)
		ndimwin = Int.(read(io, (Int32, num_kpts)))
		u_matrix_opt_data = read(io, (ComplexF64, num_bands, num_wann, num_kpts))
		inwindow_bands = Vector{Vector{Int}}(undef, num_kpts)
		u_matrix_opt = Vector{Matrix{ComplexF64}}(undef, num_kpts)
		for ik in 1:num_kpts
			inwindow_bands[ik] = include_bands[lwindow[:, ik]]
			u_matrix_opt[ik] = u_matrix_opt_data[1:ndimwin[ik], :, ik]
		end
	else
		omega_invariant = 0
		lwindow = Matrix{Bool}(undef, 0, 0)
		ndimwin = Vector{Int}(undef, 0)
		inwindow_bands = Vector{Vector{Int}}(undef, num_kpts)
		for ik in 1:num_kpts
			inwindow_bands[ik] = include_bands
		end
		u_matrix_opt = Vector{Matrix{ComplexF64}}(undef, 0)
	end

	# 8. 读取核心U矩阵
	u_matrix_data = read(io, (ComplexF64, num_wann, num_wann, num_kpts))
	u_matrix = Vector{Matrix{ComplexF64}}(undef, num_kpts)
	for ik in 1:num_kpts
		u_matrix[ik] = u_matrix_data[:, :, ik]
	end

	# 9. 读取重叠矩阵m_matrix
	m_matrix = read(io, (ComplexF64, num_wann, num_wann, nntot, num_kpts))

	# 10. 读取Wannier中心
	wannier_centres = read(io, (Float64, 3, num_wann))
	wannier_centres = map(Vec3{Float64}, eachcol(wannier_centres))

	wannier_spreads = read(io, (Float64, num_wann))

	# 返回结构化数据
	return wannier90_chk(header, num_bands,
		exclude_bands, include_bands, inwindow_bands,
		lattice, rlattice,
		mp_grid, kdirect, nntot, num_wann,
		checkpoint,
		have_disentangled,
		omega_invariant, lwindow, ndimwin, u_matrix_opt,
		u_matrix, m_matrix,
		wannier_centres, wannier_spreads,
	)
end
