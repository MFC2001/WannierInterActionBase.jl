struct KeldyshGauss2DLRCorrection <: AbstractLRCorrection
	# wannier structure.
	lattice::Lattice{Float64}
	rlattice::ReciprocalLattice{Float64}
	rldot::Mat3{Float64}
	Ω::Float64
	norb::Int
	orblocat::Vector{ReducedCoordinates{Float64}}
end
function (v::KeldyshGauss2DLRCorrection)(::Val{+}, A, k::ReducedCoordinates; atol::Real = 1e-10)
	if k ⋅ (v.rldot * k) ≤ atol
		return _KeldyshGauss2DLRCorrection_wohead!(A, v)
	else
		return _KeldyshGauss2DLRCorrection_notΓ!(A, v, k)
	end
end
function (v::KeldyshGauss2DLRCorrection)(::Val{+}, A, k::ReducedCoordinates, nkorkgrid; atol::Real = 1e-10)
	if k ⋅ (v.rldot * k) ≤ atol
		return _KeldyshGauss2DLRCorrection_wihead!(A, v, nkorkgrid)
	else
		return _KeldyshGauss2DLRCorrection_notΓ!(A, v, k)
	end
end
