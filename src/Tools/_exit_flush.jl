function _exit_flush(ios...)
	# define a global error handler that flushes all provided IO streams upon an unhandled exception.
	cleanup = function ()
		for io in ios
			flush(io)
		end
	end
	atexit(cleanup)
end
