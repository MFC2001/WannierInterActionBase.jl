"""
	_LOG_FILES_INITIALIZED

A global set that tracks the log files that have already been **initialized** (created/overwritten)
during the current process run. Ensures that:
- First write: overwrites the old file
- Subsequent writes: append to the file
- Restarting the process resets the state automatically
"""
const _LOG_FILES_INITIALIZED = Set{String}()
"""
	get_logio(logfile) -> IOStream

Open `logfile` for writing, creating the parent directory if necessary. The file is
truncated on the first call in the current process and appended to on later calls.
"""
function get_logio(logfile)
	if logfile ∈ _LOG_FILES_INITIALIZED
		return open(logfile, "a")
	else
		push!(_LOG_FILES_INITIALIZED, logfile)
		mkpath(dirname(logfile))
		return open(logfile, "w")
	end
end
