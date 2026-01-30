export ReadEIGENVAL, BANDGAP

function ReadEIGENVAL(file::AbstractString; kweight=false)
    #Scan EIGENVAL.
    file = open(file, "r")
    for _ in 1:5
        readline(file)
    end
    line = readline(file)
    line = parse.(Int, split(line))

    sumelec = line[1]
    Nk = line[2]
    Nbands = line[3]

    kpoint = Matrix{Float64}(undef, Nk, 4)
    E = Matrix{Float64}(undef, Nbands, Nk)
    occupy = Matrix{Float64}(undef, Nbands, Nk)

    T = Matrix{Float64}(undef, Nbands, 3)
    for i in 1:Nk
        kpoint[i, :] = parse.(Float64, split(myreadline(file)))

        for j in 1:Nbands
            T[j, :] = parse.(Float64, split(myreadline(file)))
        end
        E[:, i] = T[:, 2]
        occupy[:, i] = T[:, 3]
    end
    close(file)

    if kweight
        for k in 1:Nk
            if kpoint[k, 4] > eps()
                E[1, k] = NaN
            end
        end
        I = .!isnan.(E[1, :])
        E = E[:, I]
    end

    return E, kpoint
end

function myreadline(file::IOStream)
    local line
    while !eof(file)
        line = readline(file)
        if occursin(r"^\s*$", line)
            continue
        end
        break
    end
    return line
end

function BANDGAP(BANDGAPfile::AbstractString)
    #Scan BANDGAP, return fermienergy & band gap.
    file = open(BANDGAPfile, "r")
    output = [NaN, NaN]
    for i in 1:5
        line = readline(file)
        if i == 2
            output[1] = parse(Float64, split(line, ":")[2])#Band gap.
        elseif i == 5
            output[2] = parse(Float64, split(line, ":")[2])#Fermienergy.
        end
    end
    close(file)

    bandgap = output[1]
    fermienergy = output[2]

    return fermienergy, bandgap
end