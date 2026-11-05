function Base.write(io::IO, orbital::wannier90_centres, ::Type{wannier90_centres}; comment="")
    now = Dates.format(Dates.now(), "yyyy/m/d H:M:S")
    comment = comment * " " * now

    N = string(numorb(orbital) + length(orbital.atom_location))
    println(io, N)
    println(io, comment)

    if isempty(orbital.index)
        for i in eachindex(orbital.location)
            @printf(io, "%-5s\t%15.8f\t%15.8f\t%15.8f\n", "X", orbital.location[i]...)
        end
    else
        for i in eachindex(orbital.location)
            @printf(io, "%-5s\t%15.8f\t%15.8f\t%15.8f\n", "X$(orbital.index[i])", orbital.location[i]...)
        end
    end

    if !isempty(orbital.atom_location)
        if isempty(orbital.atom_name)
            for i in eachindex(orbital.atom_location)
                @printf(io, "%-5s\t%15.8f\t%15.8f\t%15.8f\n", "A", orbital.atom_location[i]...)
            end
        else
            for i in eachindex(orbital.atom_location)
                @printf(io, "%-5s\t%15.8f\t%15.8f\t%15.8f\n", "$(orbital.atom_name[i])", orbital.atom_location[i]...)
            end
        end

    end
    return nothing
end
