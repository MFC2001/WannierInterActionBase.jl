"""
	MonkhorstPack <: AbstractBrillouinZone

Perform BZ sampling employing a Monkhorst-Pack grid.

```julia
julia> kgrid = MonkhorstPack([5, 5, 1]; kshift = [1//2, 0, 0])
MonkhorstPack([5, 5, 1], [0.5, 0.0, 0.0])

julia> length(kgrid)
25

```
"""
struct MonkhorstPack <: AbstractBrillouinZone
	kgrid_size::Vec3{Int}
	kshift::Vec3{Rational{Int}}
	function MonkhorstPack(size, kshift)
		map(kshift) do ks
			ks in (0, 1 // 2) || error("Only kshifts of 0 or 1//2 implemented.")
		end
		new(size, kshift)
	end
end
MonkhorstPack(kgrid_size::AbstractVector; kshift = [0, 0, 0]) = MonkhorstPack(kgrid_size, kshift)
MonkhorstPack(kgrid_size::Tuple; kshift = [0, 0, 0]) = MonkhorstPack(kgrid_size, kshift)
MonkhorstPack(k1::Integer, k2::Integer, k3::Integer) = MonkhorstPack([k1, k2, k3])
function Base.show(io::IO, kgrid::MonkhorstPack)
	print(io, "MonkhorstPack(", kgrid.kgrid_size)
	if !iszero(kgrid.kshift)
		print(io, ", ", Float64.(kgrid.kshift))
	end
	print(io, ")")
end
Base.length(kgrid::MonkhorstPack) = prod(kgrid.kgrid_size)
