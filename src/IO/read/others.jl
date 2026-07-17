export trace, gf2dos

function trace(file::AbstractString, noperate::Integer)
    #read trace.txt from vasp2trace
    file = open(file, "r")
    nelectron = parse(Int, readline(file))
    issoc = parse(Int, readline(file))
    Noperation = parse(Int, readline(file))
    for i in 1:Noperation
        readline(file)
    end
    Nkpoint = parse(Int, readline(file))
    kpoint = zeros(Nkpoint, 3)
    for i in 1:Nkpoint
        kpoint[i, :] = parse.(Float64, split(readline(file)))
    end

    data = Vector{Pair}(undef, Nkpoint)
    aimdata = Vector{Pair}(undef, Nkpoint)
    for i in 1:Nkpoint
        noperate_ik = parse(Int, readline(file))
        operate_index = parse.(Int, split(readline(file)))
        datai = Matrix{Float64}(undef, 0, 3 + 2 * length(operate_index))
        while !eof(file)
            mark(file)
            line = readline(file)
            if occursin(r"^\s*\d*\s*$", line)
                reset(file)
                break
            else
                datai = [datai; reshape(parse.(Float64, split(line)), 1, :)]
            end
        end
        data[i] = kpoint[i, :] => datai

        I = findall(operate_index .== noperate)[1]
        aimdata[i] = kpoint[i, :] => datai[:, (2*I-1:2*I).+3]
    end
    close(file)

    Z4 = 0
    for i in Nkpoint
        Z4 += sum(aimdata[i].second[:, 1]) / 2
    end
    Z4 = Z4 % 4

    return Z4
end
#############################################################################################
function gf2dos(gf::AbstractArray{<:Real}; orbitalindex::AbstractVector{<:Pair}=["all" => 1:size(gf, 3)])
    #size(gf) = ωnum,Nk,norb
    (ωnum, Nk, _) = size(gf)

    norb = length(orbitalindex)
    dos = Vector{Pair}(undef, norb)
    for k in 1:norb
        index = orbitalindex[k].second
        T = Matrix{eltype(gf)}(undef, ωnum, Nk)
        for i in 1:Nk, j in 1:ωnum
            T[j, i] = log(sum(gf[j, i, index]))
        end
        dos[k] = orbitalindex[k].first => T
    end

    return dos
end
function greenfunction(file::AbstractString)
    #Read file(.dat) & file.

    datfile = open(file * ".dat", "r")
    readline(datfile)

    line = readline(datfile)
    line = split(line, ':')[2]
    line = parse.(Int, split(line))
    Nk = line[1]
    ωnum = line[2]

    #kindex.
    line = readline(datfile)
    line = split(line, ':')[2]
    kindex = parse.(Int, split(line))
    #kname.
    line = readline(datfile)
    line = split(line, ':')[2]
    kname = split(line)
    #greenfunction type.
    line = readline(datfile)
    gftype = match(r"^# gftype:\s*(\S+)$", line).captures[1]
    gftype = eval(Symbol(gftype))

    #kline & ω
    data = readdlm(datfile; comments=true, comment_char='#')
    close(datfile)

    kdata = KData(; line=data[1:Nk], name=kname, index=kindex)
    ω = data[Nk+1:end]

    #greenfunction
    gf = open(file, "r") do file
        return read(file)
    end

    gf = reshape(reinterpret(gftype, gf), ωnum, Nk, :)

    return kdata, ω, gf
end