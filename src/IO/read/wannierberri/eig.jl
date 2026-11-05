function Base.read(path::AbstractString, ::Type{wannierberri_eig})

    data = npzread(path)

    nk = data["NK"]
    nband = length(data["data_0"])

    band = Matrix{Float64}(undef, nband, nk)
    for ik in 1:nk
        band[:, ik] = data["data_$(ik-1)"]
    end
    return band
end
