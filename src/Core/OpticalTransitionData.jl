export OpticalTransitionData, absorption, optcond
"""
	OpticalTransitionData

Vertical optical transition data: the transition energies `ΔE` (eV), the Berry
connection matrix elements `berryconnection` (Å), the number of spin and spinor
components, the crystal volume `NΩ` (Å³) and, for 2D systems, the layer thickness
`Lz` (Å).
"""
struct OpticalTransitionData
    NΩ::Float64 # Å³
    Lz::Float64 # Å, used to calculate optical conduction of 2D. 
    nspin::Int
    nspinor::Int
    ΔE::Vector{Float64} # eV
    berryconnection::Vector{Vec3{ComplexF64}} # Å
end
function OpticalTransitionData(ΔE::AbstractVector{<:Real}, berryconnection::AbstractMatrix{<:Number};
    nspin::Integer=1, nspinor::Integer=1,
    nk::Integer=1, Ω::Real=1, NΩ::Real=nk * Ω, Lz::Real=-1)
    @assert size(berryconnection, 1) == 3 "Expected the first dimension of `berryconnection` to be 3."
    @assert length(ΔE) == size(berryconnection, 2) "Expected `length(ΔE)` to match the second dimension of `berryconnection`."
    berryconnection = map(v->Vec3{ComplexF64}(v), eachcol(berryconnection))
    return OpticalTransitionData(NΩ, Lz, nspin, nspinor, ΔE, berryconnection)
end
function OpticalTransitionData(ΔE::AbstractVector{<:Real}, berryconnection::AbstractVector{<:AbstractVector{<:Number}};
    nspin::Integer=1, nspinor::Integer=1,
    nk::Integer=1, Ω::Real=1, NΩ::Real=nk * Ω, Lz::Real=-1)
    @assert length(ΔE) == length(berryconnection) "Expected `length(ΔE)` to match the length of `berryconnection`."
    berryconnection = map(v->Vec3{ComplexF64}(v), berryconnection)
    return OpticalTransitionData(NΩ, Lz, nspin, nspinor, ΔE, berryconnection)
end
"""
	absorption(opt::OpticalTransitionData, ωgrid, ê=[0, 0, 1];
		broadening::Symbol=:gauss, fwhm=0.1, outfile=nothing) -> Vector
    absorption(opt::OpticalTransitionData, ωgrid::AbstractVector{<:Real}, μ::Integer, ν::Integer;
        broadening::Symbol=:gauss, fwhm=0.1, outfile=nothing) -> Vector

Return the imaginary part of the dielectric function ε₂(ω) on `ωgrid` (eV) from the
vertical transitions in `opt`, for the light polarization `ê` in Cartesian coordinates
or, in the second method, for the tensor component `(μ, ν)`. Transitions are broadened
with the kernel `broadening` and width `fwhm`; if `outfile` is given the result is
written to that file as well.
"""
function absorption(opt::OpticalTransitionData, ωgrid::AbstractVector{<:Real}, ê=[0, 0, 1];
    broadening::Symbol=:gauss, fwhm::Union{Real,AbstractVector{<:Real}}=0.1,
    outfile::Union{AbstractString,Nothing}=nothing)

    # ê should be Cartesian Coordinates.
    ê = ê ./ norm(ê)

    NE = length(opt.ΔE)
    coff = 1e3 * (π * qₑ) / (opt.NΩ * ϵ₀)
    if opt.nspin == opt.nspinor == 1
        coff *= 2
    end

    ϵ₂ = Vector{Float64}(undef, NE)
    Threads.@threads for i in Base.OneTo(NE)
        ϵ₂[i] = coff * abs2(ê ⋅ opt.berryconnection[i])
    end

    result = broaden_peaks(broadening, ωgrid, opt.ΔE, ϵ₂; fwhm, normalize=true)
    isnothing(outfile) || _write_absorption(outfile, ωgrid, result)

    return result
end
function absorption(opt::OpticalTransitionData, ωgrid::AbstractVector{<:Real}, μ::Integer, ν::Integer;
    broadening::Symbol=:gauss, fwhm::Union{Real,AbstractVector{<:Real}}=0.1,
    outfile::Union{AbstractString,Nothing}=nothing)

    NE = length(opt.ΔE)
    coff = 1e3 * (π * qₑ) / (opt.NΩ * ϵ₀)
    if opt.nspin == opt.nspinor == 1
        coff *= 2
    end

    ϵ₂ = Vector{Float64}(undef, NE)
    Threads.@threads for i in Base.OneTo(NE)
        ϵ₂[i] = coff * opt.berryconnection[i][μ] * conj(opt.berryconnection[i][ν])
    end

    result = broaden_peaks(broadening, ωgrid, opt.ΔE, ϵ₂; fwhm, normalize=true)
    isnothing(outfile) || _write_absorption(outfile, ωgrid, result)

    return result
end
function _write_absorption(file, ωgrid, data)
    mkpath(dirname(file))
    io = open(file, "w")
    @printf(io, "# %16s %16s\n", "ω(eV)", "ϵ₂(ω)")
    for (ω, ϵ₂) in zip(ωgrid, data)
        @printf(io, "%16.9f %16.9f\n", ω, ϵ₂)
    end
    close(io)
    return nothing
end
"""
	optcond(opt::OpticalTransitionData, ωgrid, ê=[0, 0, 1];
		broadening=:gauss, fwhm=0.1, outfile=nothing) -> Vector
    optcond(opt::OpticalTransitionData, ωgrid, μ::Integer, ν::Integer;
		broadening=:gauss, fwhm=0.1, outfile=nothing) -> Vector

Return the optical conductivity σ(ω) on `ωgrid` (eV) for the polarization `ê` or, in the
second method, for the tensor component `(μ, ν)`. The result is in S/m for 3D systems
and in units of ``e^2/\\hbar`` for 2D systems (`opt.Lz > 0`); transitions are broadened as
in [`absorption`](@ref).
"""
function optcond(opt::OpticalTransitionData, ωgrid::AbstractVector{<:Real}, ê=[0, 0, 1];
    broadening::Symbol=:gauss, fwhm::Union{Real,AbstractVector{<:Real}}=0.1,
    outfile::Union{AbstractString,Nothing}=nothing)

    # ê should be Cartesian Coordinates.
    ê = ê ./ norm(ê)

    NE = length(opt.ΔE)
    if opt.Lz > 0 # 2D
        coff = (π * opt.Lz) / opt.NΩ # e²/ħ
    else
        coff = 1e6 * (π * qₑ^2) / (opt.NΩ * ħ) # S/m
    end
    if opt.nspin == opt.nspinor == 1
        coff *= 2
    end

    σ = Vector{Float64}(undef, NE)
    Threads.@threads for i in Base.OneTo(NE)
        σ[i] = coff * abs2(ê ⋅ opt.berryconnection[i])
    end

    result = broaden_peaks(broadening, ωgrid, opt.ΔE, σ; fwhm, normalize=true)
    result .*= ωgrid

    isnothing(outfile) || begin
        if opt.Lz > 0
            _write_optcond_2D(outfile, ωgrid, result)
        else
            _write_optcond_3D(outfile, ωgrid, result)
        end
    end
    return result
end
function optcond(opt::OpticalTransitionData, ωgrid::AbstractVector{<:Real}, μ::Integer, ν::Integer;
    broadening::Symbol=:gauss, fwhm::Union{Real,AbstractVector{<:Real}}=0.1,
    outfile::Union{AbstractString,Nothing}=nothing)

    # ê should be Cartesian Coordinates.
    ê = ê ./ norm(ê)

    NE = length(opt.ΔE)
    if opt.Lz > 0 # 2D
        coff = (π * opt.Lz) / opt.NΩ # e²/ħ
    else
        coff = 1e6 * (π * qₑ^2) / (opt.NΩ * ħ) # S/m
    end
    if opt.nspin == opt.nspinor == 1
        coff *= 2
    end

    σ = Vector{Float64}(undef, NE)
    Threads.@threads for i in Base.OneTo(NE)
        σ[i] = coff * opt.berryconnection[i][μ] * conj(opt.berryconnection[i][ν])
    end

    result = broaden_peaks(broadening, ωgrid, opt.ΔE, σ; fwhm, normalize=true)
    result .*= ωgrid

    isnothing(outfile) || begin
        if opt.Lz > 0
            _write_optcond_2D(outfile, ωgrid, result)
        else
            _write_optcond_3D(outfile, ωgrid, result)
        end
    end
    return result
end
function _write_optcond_3D(file, ωgrid, data)
    mkpath(dirname(file))
    io = open(file, "w")
    @printf(io, "# %16s %16s\n", "ω(eV)", "σ(ω) (S/m)")
    for (ω, σ) in zip(ωgrid, data)
        @printf(io, "%16.9f %16.9f\n", ω, σ)
    end
    close(io)
    return nothing
end
function _write_optcond_2D(file, ωgrid, data)
    mkpath(dirname(file))
    io = open(file, "w")
    @printf(io, "# %16s %16s\n", "ω(eV)", "σ(ω) (e²/ħ)")
    for (ω, σ) in zip(ωgrid, data)
        @printf(io, "%16.9f %16.9f\n", ω, σ)
    end
    close(io)
    return nothing
end
