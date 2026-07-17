function Plots.plot(U::HR, W::HR,
	lattice::Union{Lattice, AbstractMatrix},
	orblocat::AbstractVector{<:AbstractVector{<:Real}};
	ϵ::Union{Real, Nothing} = nothing,
	title = "",
	linewidth = 1.2,
	linecolor = :black,
	ylabel = "E/eV",
	xlims = kline.line[[1, end]] + [0, 0.001],
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
