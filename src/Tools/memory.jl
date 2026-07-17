export format_memory
"""
	format_memory(bytes::Int) -> String
"""
function format_memory(bytes::Int)
	bytes ≤ 0 && return "0 B"
	units = ["B", "KB", "MB", "GB", "TB"]
	scale = min(floor(Int, log2(bytes)/10), 4)  # 最大到 TB
	value = bytes / (1024^scale)
	return @sprintf("%.2f %s", value, units[scale+1])
end
# calculate_array_memory(T::Type, dims::Tuple) = sizeof(T) * prod(dims)
