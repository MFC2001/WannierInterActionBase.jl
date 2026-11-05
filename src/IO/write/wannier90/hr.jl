function Base.write(io::IO, hr::AbstractHR{T}, ::Type{wannier90_hr}; comment="") where {T}
    now = Dates.format(Dates.now(), "yyyy/m/d H:M:S")
    comment = comment * ", " * now

    hr = HRFull(hr)
    hr = sort(hr)

    println(io, comment)

    @printf(io, "%12u\n", numorb(hr))
    @printf(io, "%12u\n", hr.nR)

    for i in 1:hr.nR
        @printf(io, "%5u", hr.ndegen[i])
        if mod(i, 15) == 0
            @printf(io, "\n")
        end
    end
    if mod(hr.nR, 15) ≠ 0
        @printf(io, "\n")
    end

    for (i, j, iR, v) in enumerate(hr)
        @printf(io, "%5u %5u %5u %7u %7u %14.8f %14.8f\n", hr.R[iR]..., i, j, reim(v)...)
    end

    return nothing
end
