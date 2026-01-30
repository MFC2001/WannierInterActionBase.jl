export WritePBAND
"""
    WritePBAND(kdata::KData, data::AbstractVector{LinearAlgebra.Eigen}, orbital::AbstractVector{<:Pair}, file::AbstractString; mode="w", comment="From MyWrite.PBAND.")
    WritePBAND(datafolder::AbstractString, orbital::AbstractVector{<:Pair}, file::AbstractString=joinpath(datafolder, "PBAND.dat"); bandindex="all", mode="w", comment="From MyWrite.PBAND.")
"""
function WritePBAND(kdata::KData, data::AbstractVector, orbital::AbstractVector{<:Pair}, file::AbstractString; mode="w", comment="From MyWrite.PBAND.")
    #orbital = ["name" => [index]].

    tablehead = "#K-Path          Energy"
    for orb in orbital
        tablehead *= @sprintf("%7s", orb.first)
    end
    tablehead *= @sprintf("%7s", "tot")

    #Processing data.
    Nk = kdata.N
    Nband = length(data[1].values)
    norb = length(orbital)

    band = Matrix{Float64}(undef, Nband, Nk)
    orbcomp = Array{Float64}(undef, Nband, Nk, norb + 1)
    for i = 1:Nk
        Et = data[i].values
        Vt = data[i].vectors
        index = sortperm(Et)
        band[:, i] = Et[index]
        Vt = Vt[:, index]
        for j in 1:Nband
            for k in 1:norb
                orbcomp[j, i, k] = sum(abs2, Vt[orbital[k].second, j])
            end
            orbcomp[j, i, norb+1] = sum(orbcomp[j, i, 1:norb])
        end
    end

    #Write data to file.
    PBAND_write(kdata, band, bandindex, orbcomp, tablehead, file; mode, comment)

    return nothing
end
function WritePBAND(datafolder::AbstractString, orbital::AbstractVector{<:Pair}, file::AbstractString=joinpath(datafolder, "PBAND.dat"); bandindex="all", mode="w", comment="From MyWrite.PBAND.")
    error("To be continued.")
    #file BAND.dat and folder wave_nk in datafolder.
    #orbital = ["name" => [index]].

    tablehead = "#K-Path          Energy"
    for orb in orbital
        tablehead *= @sprintf("%7s", orb.first)
    end
    tablehead *= @sprintf("%7s", "tot")

    (band, kdata) = MyScan.BAND(joinpath(datafolder, "BAND.dat"); readkname="Y")

    if bandindex == "all"
        Nband = size(band, 1)
        bandindex = collect(1:Nband)
    end

    #Processing data.
    Nband = length(bandindex)
    Nk = kdata.N
    norb = length(orbital)

    band = band[bandindex, :]
    orbcomp = Array{Float64}(undef, Nband, Nk, norb + 1)

    wavefolder = joinpath(datafolder, "wave_nk")
    Threads.@threads for i in 1:Nk
        wavefile = joinpath(wavefolder, "wave_$(i)")
        wave = abs2.(MyScan.WAVE(wavefile; bandindex))
        for j in 1:Nband
            for k in eachindex(orbital)
                orbcomp[j, i, k] = sum(wave[orbital[k].second, j])
            end
            orbcomp[j, i, end] = sum(orbcomp[j, i, 1:norb])
        end
    end

    #Write data to file.
    PBAND_write(kdata, band, bandindex, orbcomp, tablehead, file; mode, comment)

    return nothing
end
#################################################################################################
function PBAND_write(kdata::KData, band::AbstractMatrix{<:Real}, bandindex::AbstractVector{<:Integer}, orbcomp::AbstractArray{<:Real}, tablehead::AbstractString, file::AbstractString; mode="w", comment="From MyWrite.PBAND.")
    #Write data to PBAND file.
    now = Dates.format(Dates.now(), "yyyy/m/d H:M:S")
    comment = "#" * comment * " " * now

    (Nband, Nk) = size(band)
    norb = size(orbcomp, 3)

    path = dirname(file)
    mkpath(path)
    file = open(file, mode)

    println(file, comment)
    println(file, tablehead)
    println(file, "# NKPTS & NBANDS: ", Nk, Nband)

    write(file, "# kindex: ")
    writedlm(file, reshape(kdata.index, 1, :), "\t")
    write(file, "# kpoint: ")
    for index in kdata.index
        @printf(file, "%8.4f", kdata.line[index])
    end
    println(file, "")
    write(file, "# kname: ")
    writedlm(file, reshape(kdata.name, 1, :), "\t")


    if Nk == 1
        write(file, "# Cluster \n")
        output = Vector{String}(undef, Nband)

        Threads.@threads for i in 1:Nband
            output[i] = @sprintf("%8u%14.6f ", i, band[i])
            for j in 1:norb
                output[i] = output[i] * @sprintf("%7.3f", orbcomp[i, 1, j])
            end
        end
        writedlm(file, output)

    else
        kline = Vector{String}(undef, Nk)
        for i in 1:Nk
            kline[i] = @sprintf("%10.5f ", kdata.line[i])
        end
        output = Matrix{String}(undef, Nk, Nband)
        if Threads.nthreads() > 1
            Threads.@threads for i in 1:Nband
                for j in 1:Nk
                    output[j, i] = @sprintf("%14.6f", band[i, j])
                    for k in 1:norb
                        output[j, i] = output[j, i] * @sprintf("%7.3f", orbcomp[i, j, k])
                    end
                end
            end
            for i in 1:Nband
                println(file, "# Band-Index  ", bandindex[i])
                writedlm(file, kline .* output[:, i])
                println(file, "")
            end
        else
            for i in 1:Nband
                for j in 1:Nk
                    output[j] = @sprintf("%14.6f", band[i, j])
                    for k in 1:norb+1
                        output[j] = output[j] * @sprintf("%7.3f", orbcomp[i, j, k])
                    end
                end
                println(file, "# Band-Index  ", bandindex[i])
                writedlm(file, kline .* output)
                println(file, "")
            end
        end
    end

    close(file)
    return nothing
end
