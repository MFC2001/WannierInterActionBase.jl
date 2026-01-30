export WriteWAVE
"""
    WriteWAVE(data::LinearAlgebra.Eigen, file::AbstractString; mode="w", comment="From MyWrite.WAVE.")
    WriteWAVE(wave::AbstractArray{<:Number,3}, folder::AbstractString; mode="w", comment="From MyWrite.WAVE.")
"""
function WriteWAVE(data, file::AbstractString; mode="w", comment="From MyWrite.WAVE.")
    #Recommend file = "wave_$(nk)"
    now = Dates.format(Dates.now(), "yyyy/m/d H:M:S")
    comment = comment * " " * now

    wave = data.vectors

    path = dirname(file)
    mkpath(path)
    WAVE_write(file, wave; mode, comment)

    return nothing
end
function WriteWAVE(wave::AbstractArray{<:Number,3}, folder::AbstractString; mode="w", comment="From MyWrite.WAVE.")
    #(norb, nband, nk) = size(wave), Recommand folder = "***/wave_nk"
    now = Dates.format(Dates.now(), "yyyy/m/d H:M:S")
    comment = comment * " " * now

    folder = joinpath(folder, "wave_nk/")
    mkpath(folder)

    if Threads.nthreads() > 1
        Threads.@threads for nk in axes(wave, 3)
            file = joinpath(folder, "wave_$nk")
            WAVE_write(file, wave[:, :, nk]; mode, comment)
        end
    else
        for nk in axes(wave, 3)
            file = joinpath(folder, "wave_$nk")
            WAVE_write(file, wave[:, :, nk]; mode, comment)
        end
    end

    return nothing
end
function WAVE_write(file::AbstractString, wave::AbstractMatrix{<:Number}; mode="w", comment="From MyWrite.WAVE.")

    datfile = open(file * ".dat", mode)
    println(datfile, comment)
    println(datfile, "wavefile: ", file)
    (norb, nband) = size(wave)
    println(datfile, "wavesize: norb = ", norb, ", nband = ", nband)


    if eltype(wave) <: Complex
        println(datfile, "wavetype: ComplexF16")
        close(datfile)

        outtype = ComplexF16

    elseif eltype(wave) <: Real
        println(datfile, "wavetype: Float16")
        close(datfile)

        outtype = Float16
    end

    open(file, mode) do file
        write(file, convert.(outtype, wave))
    end

    return nothing
end