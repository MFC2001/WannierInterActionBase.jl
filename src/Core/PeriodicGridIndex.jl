export PeriodicGridIndex
"""
	PeriodicGridIndex(grid) -> PeriodicGridIndex

`grid` is a collection of `Integer`.

input: real value of grid point;
output: the index of equivelent point.
The equivelent point is in the range of 0:Ng-1, but we return the index which is in the range of 1:Ng.

If we set the size of grid as `Ng`, input `n`, return `mod(n, Ng) + 1`.

```julia
julia> pgi = PeriodicGridIndex((5, 8, 9))
PeriodicGridIndex((5, 8, 9))

julia> pgi([0,1,3])
CartesianIndex(1, 2, 4)

julia> pgi([-1,2,3])
CartesianIndex(5, 3, 4)

julia> pgi(:linear, [-1,2,3])
135

julia> (4-1)*8*5+(3-1)*5+5
135
```
"""
struct PeriodicGridIndex{N}
	grid::NTuple{N, Int}
	linearindices::LinearIndices{N, NTuple{N, Base.OneTo{Int}}}
end
function Base.show(io::IO, pgi::PeriodicGridIndex{N}) where {N}
	print(io, "PeriodicGridIndex(")
	show(io, pgi.grid)
	print(io, ")")
end
@inline _pgi_map(n::Integer, Ng::Integer) = mod(n, Ng) + 1
function (pgi::PeriodicGridIndex{N})(gridvalue) where {N}
	return CartesianIndex{N}(ntuple(d -> _pgi_map(gridvalue[d], pgi.grid[d]), Val(N)))
end
function (pgi::PeriodicGridIndex{N})(sym::Symbol, gridvalue) where {N}
	return pgi(Val(sym), gridvalue)
end
function (pgi::PeriodicGridIndex{N})(::Val{:linear}, gridvalue) where {N}
	return pgi.linearindices[pgi(gridvalue)]
end
function PeriodicGridIndex(grid)
	N = length(grid)
	grid = ntuple(i->Int(grid[i]), Val(N))
	linearindices = LinearIndices(ntuple(i -> 1:grid[i], Val(N)))
	return PeriodicGridIndex{N}(grid, linearindices)
end
