"""
	scissorslinearfit(band::AbstractMatrix{<:Real}, vindex, cindex;	cvfit::AbstractVector{<:Real} = zeros(Int, 6)) -> similar(band, Float64)

cvfit = [evs, ev0, evdel, ecs, ec0, ecdel].
ev_cor = ev_in + evs + evdel * (ev_in - ev0);
ec_cor = ec_in + ecs + ecdel * (ec_in - ec0).
You can also specify the parameters separately through `evs = 0.1`.
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
	newband[vindex, :] .= band[vindex, :] .+ evs .+ evdel .* (band[vindex, :] .- ev0)
	newband[cindex, :] .= band[cindex, :] .+ ecs .+ ecdel .* (band[cindex, :] .- ec0)

	return newband
end
