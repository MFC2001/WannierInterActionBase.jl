function parentdir(dir)
	(dir1, dir2) = splitdir(dir)
	return isempty(dir2) ? dirname(dir1) : dir1
end
