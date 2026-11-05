function Plots.plot(kline::Kline, band::AbstractVector{<:Eigen}; kwargs...)
    Nk = length(band)
    Nband = length(band[1].values)
    bandmatrix = Matrix{Float64}(undef, Nband, Nk)
    for i in eachindex(band)
        bandmatrix[:, i] .= band[i].values
    end
    return Plots.plot(kline, bandmatrix; kwargs...)
end
function Plots.plot(kline::Kline, band::AbstractMatrix{<:Real};
    weights=nothing,
    fermienergy=nothing,
    title="",
    linewidth=1.2,
    linecolor=:black,
    ylabel="E/eV",
    xlims=kline.line[[1, end]] + [-0.001, 0.001],
    ylims=(minimum(band), maximum(band)),
    xticks=(kline.line[kline.index], kline.name),
    legend=false,
    size=(600, 700),
	scale = 1,
	markercolor=:red,
    alpha=0.5,
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
    if !isnothing(weights)
        weights = _bandplot_format_weights(weights)
        (nband, nk) = Base.size(weights)
        @assert Base.size(band) == (nband, nk) "Wrong size of weights!"
        # 1D data
        xs = repeat(kline.line, inner=nband)
        ys = vec(band)
        ms = vec(weights) .* scale .+ 1e-6 # zero will give a dot with finite size
        Plots.scatter!(p, xs, ys;
            markersize=ms,
            markercolor,
            markerstrokewidth=0,
            alpha,
            label=nothing
        )
    end
    if !isnothing(fermienergy)
        Plots.plot!(kline.line[[1, end]], fermienergy * ones(2); linecolor=:black, linestyle=:dot)
    end
    return p
end
function _bandplot_format_weights(weights::AbstractVector{<:AbstractVector{<:Real}})
    nb = length(weights[1])
    nk = length(weights)
    weights_mat = Matrix{eltype(weights[1])}(undef, nb, nk)
    for ik in 1:nk
        weights_mat[:, ik] .= weights[ik]
    end
    return weights_mat
end
function _bandplot_format_weights(weights::AbstractVector{<:Real})
    nk = length(weights)
    weights_mat = Matrix{eltype(weights[1])}(undef, 1, nk)
    for ik in 1:nk
        weights_mat[1, ik] = weights[ik]
    end
    return weights_mat
end
function _bandplot_format_weights(weights::AbstractMatrix{<:Real})
    return weights
end
