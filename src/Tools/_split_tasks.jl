"""
	_split_tasks(all_tasks, nchunks::Integer) -> Vector{UnitRange{Int}}
	_split_tasks(ntasks::Integer, nchunks::Integer) -> Vector{UnitRange{Int}}
"""
# 虽然Julia的多线程会自动进行分块，但是如果有一些缓存数据需要预分配内存并绑定到任务块以减少内存分配，手动进行任务分块是有必要的。
function _split_tasks(all_tasks, nchunks::Integer)
	return map(idx->all_tasks[idx], _split_tasks(length(all_tasks), nchunks))
end
function _split_tasks(ntasks::Integer, nchunks::Integer)
	ntasks = Int(ntasks)
	nchunks = Int(nchunks)

	@assert ntasks > 0 "Number of tasks must be positive."
	@assert nchunks > 0 "Number of chunks must be positive."

	chunk_size = ntasks ÷ nchunks
	remainder = ntasks % nchunks

	task_chunks = Vector{UnitRange{Int}}(undef, nchunks)
	start_idx = 1
	for i in 1:nchunks
		# 前remainder个分片多1个任务，保证均匀
		end_idx = start_idx - 1 + chunk_size + (i ≤ remainder ? 1 : 0)
		task_chunks[i] = start_idx:end_idx
		start_idx = end_idx + 1
	end
	return task_chunks
end
