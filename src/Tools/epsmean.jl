function epsmean(::Type{epsmat_yambo}, outfile::AbstractString, epsfile::AbstractString, epsfiles::AbstractString ...;
	format::Symbol = :netcdf4,
)
	eps = NCDataset(epsfile; format)
	iq = _yambo_eps_get_iq(keys(eps))
	FREQ_PARS_sec_iq = eps["FREQ_PARS_sec_iq$iq"]
	FREQ_PARS_sec_iq_dims = dimnames(FREQ_PARS_sec_iq)
	FREQ_PARS_sec_iq = FREQ_PARS_sec_iq[:]
	FREQ_sec_iq = eps["FREQ_sec_iq$iq"]
	FREQ_sec_iq_dims = dimnames(FREQ_sec_iq)
	FREQ_sec_iq = FREQ_sec_iq[:, :]
	epsmat = eps["X_Q_$iq"]
	eps_dims = dimnames(epsmat)
	epsmat = epsmat[:, :, :, :]
	close(eps)
	for epsfile in epsfiles
		eps_ = NCDataset(epsfile; format)
		iq_ = _yambo_eps_get_iq(keys(eps_))
		epsmat += eps_["X_Q_$iq_"][:, :, :, :]
		close(eps_)
	end

	neps = length(epsfiles) + 1
	epsmat ./= neps

	mkpath(dirname(outfile))
	out = NCDataset(outfile, "c"; format)
	defVar(out, "FREQ_PARS_sec_iq$iq", FREQ_PARS_sec_iq, FREQ_PARS_sec_iq_dims)
	defVar(out, "FREQ_sec_iq$iq", FREQ_sec_iq, FREQ_sec_iq_dims)
	defVar(out, "X_Q_$iq", epsmat, eps_dims)
	close(out)

	return nothing
end
function _yambo_eps_get_iq(varnames)
	pattern = r"^FREQ_PARS_sec_iq(\d+)$"
	iq = nothing
	for varname in varnames
		m = match(pattern, varname)
		if !isnothing(m)
			iq = parse(Int, m.captures[1])
			break
		end
	end
	return iq
end
function epsmean(::Type{epsmat_bgw}, outfile::AbstractString, epsfile::AbstractString, epsfiles::AbstractString ...;)

end
