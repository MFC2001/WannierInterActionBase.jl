function Plots.plot(int::PairIntHR, lattice; kwargs...)
    lattice = parent(lattice)
    ldot = Mat3{Float64}(transpose(lattice) * lattice)
    normR = Vector{Float64}(undef, int.nR)
    for i in eachindex(int.R)
        normR[i] = sqrt(int.R[i] ⋅ (ldot * int.R[i]))
    end
    max_dict = Dict{Int,Float64}()
    for (ir, va) in zip(int.iRs, int.vals)
        abs_va = abs(va)
        if haskey(max_dict, ir)
            if abs_va > max_dict[ir]
                max_dict[ir] = abs_va
            end
        else
            max_dict[ir] = abs_va
        end
    end
    r_points = [normR[ir] for ir in keys(max_dict)]
    v_points = [max_dict[ir] for ir in keys(max_dict)]
    p = Plots.plot(r_points, v_points;
        seriestype=:scatter,
        kwargs...,
    )
    return p
end
