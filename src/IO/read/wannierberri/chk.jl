function Base.read(path::AbstractString, ::Type{wannierberri_chk}; atol::Real = sqrt(eps(Float64)))

	data = npzread(path)

	# 获取当前时间（UTC）
	now_time = Dates.now()
	# 拼接日期部分：日（两位） + 月缩写 + 四位年份
	date_str = Dates.format(now_time, "dd") * Dates.format(now_time, "u") * Dates.format(now_time, "yyyy")
	# 时间部分：时:分:秒（24小时制，两位）
	time_str = Dates.format(now_time, "HH:MM:SS")
	# 输出完整字符串
	header = "read on $date_str at $time_str"

	# nband
	num_bands = data["num_bands"]
	exclude_bands = Int[]
	include_bands = collect(1:num_bands)

	# lattice
	lattice = Lattice(transpose(data["real_lattice"]))
	rlattice = reciprocal(lattice)

	# kgrid
	num_kpts = data["num_kpts"]
	mp_grid = data["mp_grid"]
	kdirect = data["kpt_red"]
	kdirect = map(ReducedCoordinates{Float64}, eachrow(kdirect))

	# others
	nntot = 0
	num_wann = data["num_wann"]
	checkpoint = "postwann"
	have_disentangled = true

	omega_invariant = 0

	# 核心U矩阵
	lwindow = falses(num_bands, num_kpts)
	ndimwin = Vector{Int}(undef, num_kpts)
	inwindow_bands = Vector{Vector{Int}}(undef, num_kpts)
	u_matrix_total = Vector{Matrix{ComplexF64}}(undef, num_kpts)
	for ik in 1:num_kpts
		u_matrix_ik = data["v_matrix_$(ik-1)"]
		inwindow_bands[ik] = findall(x -> any(xx->abs(xx)>atol, x), eachrow(u_matrix_ik))
		ndimwin[ik] = length(inwindow_bands[ik])
		lwindow[inwindow_bands[ik], ik] .= true
		u_matrix_total[ik] = u_matrix_ik[inwindow_bands[ik], :]
	end
	u_matrix_opt = u_matrix_total
	
	u_matrix = Vector{Matrix{ComplexF64}}(undef, num_kpts)
	unit_mat = Matrix{ComplexF64}(I, num_wann, num_wann)
	for ik in 1:num_kpts
		u_matrix[ik] = unit_mat
	end

	m_matrix = Array{ComplexF64}(undef, num_wann, num_wann, nntot, num_kpts)
	
	# Wannier中心
	wannier_centres = data["wannier_centers_cart"]
	wannier_centres = map(CartesianCoordinates{Float64}, eachrow(wannier_centres))
	wannier_centres = map(x -> lattice \ x, wannier_centres)

	wannier_spreads = data["wannier_spreads"]

	# 返回结构化数据
	return wannier90_chk(header, num_bands,
		exclude_bands, include_bands, inwindow_bands,
		lattice, rlattice,
		mp_grid, kdirect, nntot, num_wann,
		checkpoint,
		have_disentangled,
		omega_invariant, lwindow, ndimwin, 
		u_matrix_opt, u_matrix, u_matrix_total, m_matrix,
		wannier_centres, wannier_spreads,
	)
end
