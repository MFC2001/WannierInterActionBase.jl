"""
    _fold_by_grid(R, grid)

The interval is open on the left and closed on the right, e.g.(-grid/2, grid/2].
"""
function _fold_by_grid(R, grid)
    length(R) == length(grid) ||
        throw(DimensionMismatch("R and grid must have the same length"))
    return map(_fold_by_grid_one, R, grid)
end
function _fold_by_grid!(R, grid)
    length(R) == length(grid) ||
        throw(DimensionMismatch("R and grid must have the same length"))
    for (i, r) in enumerate(R)
        R[i] = _fold_by_grid_one(r, grid[i])
    end
    return R
end
@inline function _fold_by_grid_one(x::Integer, g::Integer)
    g > 0 || throw(ArgumentError("grid element must be positive"))
    r = mod(x, g)
    return 2 * r > g ? r - g : r
end
