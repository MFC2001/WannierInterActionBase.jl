export Writelammps

function Writelammps(data, file::AbstractString, mode="w"; comment="From MyWrite.lmpdata.")
    now = Dates.format(Dates.now(), "yyyy/m/d H:M:S")
    comment = comment * " " * now

    path = dirname(file)
    mkpath(path)
    file = open(file, mode)
    write(file, comment * "\n")
    writeheader(file, data)
    writebody(file, data)
    close(file)
end
function writeheader(file::IOStream, data)
    write(file, "#Header\n")
    header = headerdict()
    for key in header["outorder"]
        if !haskey(data, key)
            continue
        end
        write(file, "\t")
        for i in eachindex(data[key])
            write(file, "$(data[key][i])\t")
        end
        write(file, key, "\n")
    end
    return nothing
end
function writebody(file::IOStream, data)
    write(file, "#Body\n")
    body = bodydict()
    for key in body["outorder"]
        if haskey(data, key) && typeof(data[key]) ≠ Symbol
            eval(body[key])(file, data)
        end
    end
    return nothing
end
function Atoms(file::IOStream, data)
    write(file, "Atoms")
    if haskey(data, "atom_style")
        write(file, " # $(data["atom_style"])\n\n")
    else
        write(file, "\n\n")
    end
    atoms = data["Atoms"]
    for i in eachindex(atoms.index)
        write(file, "\t$(atoms.index[i])\t$(atoms.type[i])")
        for j in 1:3
            write(file, "\t$(atoms.location[i,j])")
        end
        write(file, "\n")
    end
    write(file, "\n")
end
function Masses(file::IOStream, data)
    write(file, "Masses\n\n")
    masses = data["Masses"]
    for i in eachindex(masses.type)
        write(file, "\t$(masses.type[i])\t$(masses.value[i])\n")
    end
    write(file, "\n")
end