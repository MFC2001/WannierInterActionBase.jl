"""
    scissorslinearfit(band::AbstractMatrix{<:Real}, vindex, cindex;
        cvfit::AbstractVector{<:Real} = zeros(6)) -> similar(band, Float64)

Apply a scissors operator with a linear energy-dependent correction to `band`.

The correction is
- `ev_cor = ev_in + evs + evdel * (ev_in - ev0)` for valence bands,
- `ec_cor = ec_in + ecs + ecdel * (ec_in - ec0)` for conduction bands.

The six parameters are given by `cvfit = [evs, ev0, evdel, ecs, ec0, ecdel]`. 
They may also be specified individually via the keyword arguments `evs, ev0, evdel, ecs, ecdel`, which override the corresponding entries of `cvfit`.
"""
function scissorslinearfit(band::AbstractMatrix{<:Real}, vindex, cindex;
	cvfit::AbstractVector{<:Real} = zeros(Int, 6),
	evs::Real = cvfit[1],
	ev0::Real = cvfit[2],
	evdel::Real = cvfit[3],
	ecs::Real = cvfit[4],
	ec0::Real = cvfit[5],
	ecdel::Real = cvfit[6],
)

	newband = similar(band, Float64)
	isempty(vindex) || (newband[vindex, :] .= band[vindex, :] .+ evs .+ evdel .* (band[vindex, :] .- ev0))
	isempty(cindex) || (newband[cindex, :] .= band[cindex, :] .+ ecs .+ ecdel .* (band[cindex, :] .- ec0))

	return newband
end
