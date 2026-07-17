export greenfunction, gf2dos, dos

function greenfunction(kdata::KData, ω::AbstractVector{<:Real}, greenfunctiondata::AbstractArray{<:Real}, file::AbstractString; mode="w", comment="From MyWrite.dos.")
    #size(greenfunctiondata) = (ωnum,Nk,norb)
    now = Dates.format(Dates.now(), "yyyy/m/d H:M:S")
    comment = "#" * comment * " " * now

    kstring = similar(kdata.line, String)
    for i in 1:kdata.N
        kstring[i] = @sprintf("%15.6E", kdata.line[i])
    end

    ωstring = similar(ω, String)
    for i in eachindex(ω)
        ωstring[i] = @sprintf("%15.6E", ω[i])
    end

    path = dirname(file)
    mkpath(path)
    datfile = open(file * ".dat", mode)

    println(datfile, comment)
    println(datfile, "# NKPTS & ωnum: ", kdata.N, " ", length(ω))
    write(datfile, "# kindex: ")
    writedlm(datfile, reshape(kdata.index, 1, :), "\t")
    write(datfile, "# kname: ")
    writedlm(datfile, reshape(kdata.name, 1, :), "\t")

    println(datfile, "# gftype: Float32")

    println(datfile, "#    k(1/Ang)")
    writedlm(datfile, kstring)
    println(datfile, "#     E(eV)")
    writedlm(datfile, ωstring)

    close(datfile)

    open(file, mode) do file
        write(file, convert.(Float32, greenfunctiondata))
    end

    return nothing
end
#####################################################################################################
function gf2dos(gffile::AbstractString, dosfile::AbstractString; orbitalindex="all")

    (kdata, ω, gf) = MyScan.greenfunction(gffile)

    if orbitalindex == "all"
        orbitalindex = ["all" => 1:size(gf, 3)]
    elseif typeof(orbitalindex) <: AbstractVector{<:Pair}
    else
        error("Wrong orbitalindexfrom gf2dos.")
    end

    ωnum = length(ω)
    Nk = kdata.N

    norb = length(orbitalindex)
    dosdata = Vector{Pair}(undef, norb)
    for k in 1:norb
        index = orbitalindex[k].second
        T = Matrix{eltype(gf)}(undef, ωnum, Nk)
        for i in 1:Nk, j in 1:ωnum
            T[j, i] = log(sum(gf[j, i, index]))
        end
        dosdata[k] = orbitalindex[k].first => T
    end

    dos(kdata, ω, dosdata, dosfile)

    return nothing
end
###################################################################################################
function dos(kdata::KData, ω::AbstractVector{<:Real}, dosdata::AbstractVector{<:Pair}, file::AbstractString; mode="w", comment="From MyWrite.dos.")
    #dosdata = ["name" => dos].
    now = Dates.format(Dates.now(), "yyyy/m/d H:M:S")
    comment = "#" * comment * " " * now

    tablehead = "#    k(1/Ang)          E(eV)"
    for d in dosdata
        tablehead *= @sprintf("%15s", d.first)
    end

    ωstring = similar(ω, String)
    for i in eachindex(ω)
        ωstring[i] = @sprintf("%15.6E", ω[i])
    end

    Nk = kdata.N
    Nω = length(ω)
    nd = length(dosdata)

    path = dirname(file)
    mkpath(path)
    file = open(file, mode)

    println(file, comment)
    write(file, "# kindex: ")
    writedlm(file, reshape(kdata.index, 1, :), "\t")
    write(file, "# kpoint: ")
    for index in kdata.index
        @printf(file, "%8.4f", kdata.line[index])
    end
    println(file, "")
    write(file, "# kname: ")
    writedlm(file, reshape(kdata.name, 1, :), "\t")

    println(file, tablehead)

    if Nk == 1
        error("To be continued!")
    else
        if Threads.nthreads() > 1
            kstring = Vector{String}(undef, Nk)
            dosstring = fill("", Nω, Nk)
            Threads.@threads for i in 1:Nk
                kstring[i] = @sprintf("%15.6E", kdata.line[i])
                for j in 1:Nω, k in 1:nd
                    dosstring[j, i] *= @sprintf("%15.6E", dosdata[k].second[j, i])
                end
            end
            for i in 1:Nk
                writedlm(file, repeat([kstring[i]], Nω) .* ωstring .* dosstring[:, i])
                println(file, "")
            end
        else
            dosstring = Vector{String}(undef, Nω)
            for i in 1:Nk
                fill!(dosstring, "")
                kstring = @sprintf("%15.6E", kdata.line[i])
                for j in 1:Nω, k in 1:nd
                    dosstring[j] *= @sprintf("%15.6E", dosdata[k].second[j, i])
                end
                writedlm(file, repeat([kstring], Nω) .* ωstring .* dosstring)
                println(file, "")
            end
        end
    end

    close(file)
    return nothing
end