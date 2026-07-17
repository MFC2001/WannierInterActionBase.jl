export Readlammps

function Readlammps(datafile::AbstractString)
    #Read lammps data file, to be continued! Can check the dicts below.
    file = open(datafile, "r")
    readline(file) #skip the first line.
    header = readheader(file)
    body = readbody(file, header)
    close(file)
    return merge(header, body)
end

function readheader(file::IOStream)::Dict{String,Any}
    header = headerdict()
    while !eof(file)
        mark(file)
        line = readline(file)
        if occursin(r"^\s*(#|$)", line)
            continue
        end
        line = split(discard_comment(line))
        index = occursin.(r"^([a-z]|[A-Z])+$", line)
        keyword = join(line[index], " ")
        if !haskey(header, keyword)
            reset(file)
            break
        end
        index = findall(.!index)
        value = Int[]
        for i in index
            if occursin(r"^-?\d+$", line[i])
                value = [value; parse(Int, line[i])]
            elseif occursin(r"^-?\d+\.\d*$", line[i])
                value = [value; parse(Float64, line[i])]
            end
        end
        if length(value) == 1
            value = value[1]
        end
        header[keyword] = value
    end
    lohi2lattvec!(header)
    return header
end
function lohi2lattvec!(header::Dict{String,Any})
    lattvec = zeros(3, 3)
    lattvec[1, 1] = header["xlo xhi"][2] - header["xlo xhi"][1]
    lattvec[2, 2] = header["ylo yhi"][2] - header["ylo yhi"][1]
    lattvec[3, 3] = header["zlo zhi"][2] - header["zlo zhi"][1]
    lattvec[2, 1] = header["xy xz yz"][1]
    lattvec[3, 1] = header["xy xz yz"][2]
    lattvec[3, 2] = header["xy xz yz"][3]
    header["lattvec"] = lattvec
    return nothing
end

function readbody(file::IOStream, header::Dict{String,Any})::Dict{String,Any}
    body = bodydict()
    while !eof(file)
        mark(file)
        line = readline(file)
        if occursin(r"^\s*(#|$)", line)
            continue
        end
        reset(file)
        keyword = discard_comment(line, "#")
        body[keyword] = eval(body[keyword])(file)
    end
    body["atom_style"] = body["Atoms"][2]
    body["Atoms"] = body["Atoms"][1]
    (_, n) = size(body["Atoms"].location)
    if n == 6
        atomslocat = zeros(header["atoms"], 3)
        for i in 1:header["atoms"]
            atomslocat[i, :] = (body["Atoms"].location[i, 1:3] +
                                body["Atoms"].location[i, 4] * header["lattvec"][1, :] +
                                body["Atoms"].location[i, 5] * header["lattvec"][2, :] +
                                body["Atoms"].location[i, 6] * header["lattvec"][3, :])
        end
        body["Atoms"] = AtomsData(body["Atoms"].index, body["Atoms"].type, atomslocat)
    end
    return body
end

function Masses(file::IOStream)::MassesData
    #The pointer needs to be at the beginning of the line where "Masses" is located.
    line = readline(file)
    if discard_comment(line, "#") ≠ "Masses"
        error("Wrong Masses keyword!")
    end
    readline(file) #skip the line behind keyword.

    bodykey = keys(bodydict())
    massestype = Vector{Int}(undef, 0)
    massesvalue = Vector{Float64}(undef, 0)
    while !eof(file)
        mark(file)
        line = readline(file)
        if occursin(r"^\s*(#|$)", line)
            continue
        end
        line = discard_comment(line, "#")
        if line ∈ bodykey
            reset(file)
            break
        end
        line = split(line)
        push!(massestype, parse(Int, line[1]))
        push!(massesvalue, parse(Float64, line[2]))
    end
    return MassesData(massestype, massesvalue)
end

function Atoms(file::IOStream)
    #The pointer needs to be at the beginning of the line where "Atoms" is located.
    line = readline(file)
    line = match(r"^\s*(Atoms)[\s#]*([a-z]*)", line)
    if line.captures[1] ≠ "Atoms"
        error("Wrong Atoms keyword!")
    end
    atom_style = line.captures[2]
    Atom_styles = keys(atom_styledict())
    while true
        if atom_style ∈ Atom_styles
            break
        end
        print("Wrong atom_style! Please input atom_style:\n")
        atom_style = readline(stdin)
    end
    readline(file) #skip the line behind keyword.

    #judge if unitcell index exits or not.
    local n
    while !eof(file)
        mark(file)
        line = readline(file)
        if occursin(r"^\s*(#|$)", line)
            continue
        end
        line = discard_comment(line, "#")
        line = split(line)
        n = length(line) - 2
        reset(file)
        break
        if line ∈ bodykey
            reset(file)
            break
        end
    end

    bodykey = keys(bodydict())
    atomsindex = Vector{Int}(undef, 0)
    atomstype = Vector{Int}(undef, 0)
    atomslocat = Matrix{Float64}(undef, 0, n)
    while !eof(file)
        mark(file)
        line = readline(file)
        if occursin(r"^\s*(#|$)", line)
            continue
        end
        line = discard_comment(line, "#")
        if line ∈ bodykey
            reset(file)
            break
        end
        line = split(line)
        push!(atomsindex, parse(Int, line[1]))
        push!(atomstype, parse(Int, line[2]))
        atomslocat = [atomslocat; transpose(parse.(Float64, line[3:end]))]
    end
    return AtomsData(atomsindex, atomstype, atomslocat), atom_style
end

function Velocities(file::IOStream)
    #To be continued!!!
    #The pointer needs to be at the beginning of the line where "Velocities" is located.
    line = readline(file)
    if discard_comment(line, "#") ≠ "Velocities"
        error("Wrong Velocities keyword!")
    end
    readline(file) #skip the line behind keyword.

    bodykey = keys(bodydict())
    while !eof(file)
        mark(file)
        line = readline(file)
        if occursin(r"^\s*(#|$)", line)
            continue
        end
        line = discard_comment(line, "#")
        if line ∈ bodykey
            reset(file)
            break
        end
    end
    return :(Velocities)
end

function discard_comment(string::AbstractString, tag="#")
    #Discard character after tag.
    reg = Regex("^\\s*([^$tag\$]*)")
    string = match(reg, string).captures[1]
    string = replace(string, r"(\s*)$" => "")
    return string
end