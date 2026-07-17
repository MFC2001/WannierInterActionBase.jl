function Plots.plot(kline::Kline, band::AbstractVector{<:Eigen}; kwards...)
	Nk = length(band)
	Nband = length(band[1].values)
	bandmatrix = Matrix{Float64}(undef, Nband, Nk)
	for i in eachindex(band)
		bandmatrix[:, i] .= band[i].values
	end
	return Plots.plot(kline, bandmatrix; kwards...)
end

function Plots.plot(kline::Kline, band::AbstractMatrix{<:Real};
	title = "",
	linewidth = 1.2,
	linecolor = :black,
	ylabel = "E/eV",
	xlims = kline.line[[1, end]] + [-0.001, 0.001],
	ylims = (minimum(band), maximum(band)),
	xticks = (kline.line[kline.index], kline.name),
	legend = false,
	size = (600, 700),
	fermienergy = nothing,
)
	p = Plots.plot(kline.line, transpose(band);
		title,
		linewidth,
		linecolor,
		ylabel,
		xlims,
		ylims,
		xticks,
		legend,
		size,
	)
	if !isnothing(fermienergy)
		Plots.plot!(kline.line[[1, end]], fermienergy * ones(2); linecolor = :black, linestyle = :dot)
	end
	return p
end
