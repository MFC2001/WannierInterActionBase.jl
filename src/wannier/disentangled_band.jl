# 不能保证opt作用后的哈密顿量是对角的，只能保证两次幺正变换后得到的是wannier基|ik>下的哈密顿矩阵。
function disentangled_band(band::AbstractVector{<:Eigen}, u_matrix_opt::Array{<:Number,3},
    inwindow_bands::AbstractVector{<:AbstractVector{<:Integer}})

    nband = length(band[1].values)
    bandindex = setdiff(collect(1:nband), exclude_bands)

    band_opt = similar(band)

    if isempty(u_matrix_opt)
        Threads.@threads for ik in eachindex(band)
            values = band[k].values[inwindow_bands[ik]]
            vectors = band[k].vectors[:, inwindow_bands[ik]]
            band_opt[k] = Eigen(values, vectors)
        end
    else
        nwann = size(u_matrix_opt, 2)
        nband = size(u_matrix_opt, 1)
        Threads.@threads for k in eachindex(band)
            values = Vector{Float64}(undef, nwann)
            for i in 1:nwann
                values[i] = sum(n -> abs2(u_matrix_opt[n, i, k]) * band[k].values[bandindex[n]], Base.OneTo(nband))
            end
            vectors = band[k].vectors[:, bandindex] * u_matrix_opt[:, :, k]
            band_opt[k] = Eigen(values, vectors)
        end
    end

    return band_opt
end

function disentangled_band(band::AbstractVector{<:Eigen}, u_matrix_opt::Array{<:Number,3})

    if isempty(u_matrix_opt)
        band_opt = deepcopy(band)
    else
        band_opt = similar(band)

        nwann = size(u_matrix_opt, 2)
        nband = size(u_matrix_opt, 1)
        Threads.@threads for k in eachindex(band)
            values = Vector{Float64}(undef, nwann)
            for i in 1:nwann
                values[i] = sum(n -> abs2(u_matrix_opt[n, i, k]) * band[k].values[n], Base.OneTo(nband))
            end
            vectors = band[k].vectors * u_matrix_opt[:, :, k]
            band_opt[k] = Eigen(values, vectors)
        end
    end

    return band_opt
end
