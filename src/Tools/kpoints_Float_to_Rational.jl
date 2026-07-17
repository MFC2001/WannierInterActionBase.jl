function kpoints_Float_to_Rational(kpoints::Vector{<:ReducedCoordinates}, kgrid; atol::Real = 1e-6)::Vector{ReducedCoordinates{Rational{Int}}}
	k_denominator = map(1:3) do i
		kgrid.kshift[i] == 0 ? kgrid.kgrid_size[i] : kgrid.kgrid_size[i] * 2
	end
	kpoints_frac = Vector{ReducedCoordinates{Rational{Int}}}(undef, length(kpoints))
	for (ik, k) in enumerate(kpoints)
		k_numerator = k .* k_denominator
		k_numerator_Int = round.(Int, k_numerator)
		any(x->abs(x) > atol, k_numerator - k_numerator_Int) && error("Get a kpoint which is not a Rational, or get a wrong kgrid!")
		kpoints_frac[ik] = ReducedCoordinates{Rational{Int}}(k_numerator_Int .// k_denominator)
	end
	return kpoints_frac
end
function kpoints_Float_to_Rational(kpoints::Vector{<:ReducedCoordinates}; tol::Real = 1e-6)::Vector{ReducedCoordinates{Rational{Int}}}
	return map(k -> rationalize.(k; tol), kpoints)
end
function kpoints_Float_to_Rational(kpoint::ReducedCoordinates; tol::Real = 1e-6)::ReducedCoordinates{Rational{Int}}
	return rationalize.(kpoint; tol)
end
