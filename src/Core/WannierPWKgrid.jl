export WannierPWKgrid, translate!
"""
    WannierPWKgrid{WT<:Number}

Stores a set of Wannier functions expanded in a plane-wave basis.

This type is immutable: its fields cannot be changed, except for the `fft_grid` field.

See also [`QE_wannier`](@ref), [`translate!`](@ref).
"""
struct WannierPWKgrid{WT<:Number} # PWKgridStates{WT <: Number}
    nk::Int
    nw::Int
    nspinor::Int
    nG::Vector{Int} # length.(Gidx)
    kgrid::RedKgrid
    Gidx::Vector{Vector{Int}} #[nk][nG]
    Glist::Vector{ReducedCoordinates{Int}}
    components::Vector{Matrix{WT}} # [nk][nG, nw] or [nk][nG*2, nw]
    fft_grid::Ref{NTuple{3,Int}}
    Gmap::Vector{Int}
    period::Vec3{Bool}
    lattice::Lattice{Float64}
    rlattice::ReciprocalLattice{Float64}
    ldot::Mat3{Float64}
    rldot::Mat3{Float64}
    Ω::Float64
    rΩ::Float64
end
function Base.getproperty(wave::WannierPWKgrid{WT}, name::Symbol) where {WT}
    if name === :fft_grid
        # 解引用 Ref，直接返回 NTuple
        return getfield(wave, :fft_grid)[]
    elseif name === :nG_max
        return maximum(getfield(wave, :nG))
    else
        return getfield(wave, name)
    end
end
function Base.setproperty!(wave::WannierPWKgrid{WT}, name::Symbol, value) where {WT}
    if name === :fft_grid
        value = (value[1], value[2], value[3])
        getfield(wave, :fft_grid)[] = value
        pgi = PeriodicGridIndex(value)
        wave.Gmap .= pgi.(Val(:linear), wave.Glist)
        return value
    else
        error("Only allow to change the fft_grid of WannierPWKgrid!")
    end
end
"""
    WannierPWKgrid(kgrid::RedKgrid, components::Vector{Matrix{WT}},
        Gidx::Vector{Vector{Int}}, Glist::Vector{ReducedCoordinates{Int}}, lattice::Lattice;
        nspinor::Integer=1, fft_grid=(0, 0, 0), period=[true, true, true]) -> WannierPWKgrid{WT}

Construct a `WannierPWKgrid` representing Wannier functions in a plane-wave basis.

# Arguments
- `kgrid`: Full k-point grid containing all k points.
- `components`: Length-`nk` vector, where `nk` is the number of k points. 
  `components[ik]` is an `[nG_k, nw]` matrix, with `nw` the number of Wannier functions and `nG_k` the number of plane waves at k point `ik`.
- `Gidx`: Length-`nk` vector. `Gidx[ik]` is an integer vector of length `nG_k` giving indices into `Glist`. 
  Its order must match the plane-wave basis order in `components[ik]`.
- `Glist`: Complete G-vector list, i.e. the union of the plane-wave bases over all k points.
- `lattice`: Lattice information.
- `nspinor`: `1` for scalar or `2` for spinor. 
  For `nspinor == 2`, `components[ik]` has size `[2nG_k, nw]`; rows `1:nG_k` and `nG_k+1:2nG_k` correspond to the two spinor components.
- `fft_grid`: FFT grid size used in subsequent FFT-related operations.
- `period`: Periodicity along each direction. Set according to the physical system. 
  For example, for a monolayer nonperiodic in z, use `period = [true, true, false]`, even if the DFT calculation uses vacuum to maintain periodic boundary conditions.
"""
function WannierPWKgrid(kgrid::RedKgrid, components::Vector{Matrix{WT}},
    Gidx::Vector{Vector{Int}}, Glist::Vector{ReducedCoordinates{Int}}, lattice::Lattice;
    nspinor::Integer=1, fft_grid=(0, 0, 0), period=[true, true, true],
) where {WT<:Number}
    @assert (nspinor == 1 || nspinor == 2) "nspinor must be 1 or 2, got $nspinor"
    check_lattice_period(lattice, period)
    # 对于 nspinor=2 的旋量波函数，要求前半和后半为不同自旋。
    nk = length(kgrid)
    nw = size(components[1], 2)
    nG = length.(Gidx)
    fft_grid = (fft_grid[1], fft_grid[2], fft_grid[3])
    if all(iszero, fft_grid)
        Gmap = Vector{Int}(undef, length(Glist))
    else
        pgi = PeriodicGridIndex(fft_grid)
        Gmap = pgi.(Val(:linear), Glist)
    end
    period = Vec3{Bool}(collect(Bool, period))
    # lattice
    lattice = convert(Lattice{Float64}, lattice)
    rlattice = reciprocal(lattice)
    ldot = parent(lattice)
    ldot = transpose(ldot) * ldot
    rldot = parent(rlattice)
    rldot = transpose(rldot) * rldot
    Ω = abs(det(parent(lattice)))
    rΩ = abs(det(parent(rlattice)))
    return WannierPWKgrid{WT}(
        nk, nw, nspinor, nG, kgrid, Gidx, Glist, components, fft_grid, Gmap,
        period, lattice, rlattice, ldot, rldot, Ω, rΩ,
    )
end
"""
    (wave::WannierPWKgrid)(r::AbstractVector{<:Real}, R::AbstractVector{<:Integer}=ReducedCoordinates(0, 0, 0))

Evaluate the Wannier functions stored in `wave` at position `r`, for the Wannier functions centered in cell `R`.

Returns a vector of length `nw` if `wave.nspinor == 1`, or a `[2, nw]` matrix if `wave.nspinor == 2`.
"""
function (wave::WannierPWKgrid{WT})(r::AbstractVector{<:Real},
    R::AbstractVector{<:Integer}=ReducedCoordinates(0, 0, 0)) where {WT}
    r = ReducedCoordinates(r)
    eGr = map(G->cispi(2 * (G ⋅ r)), wave.Glist)
    ekr = map(k->cispi(2 * (k ⋅ r)), wave.kgrid) ./ sqrt(wave.Ω)
    ekR = map(k->cispi(2 * (k ⋅ R)), wave.kgrid) ./ wave.nk
    if wave.nspinor == 1
        ψr = Vector{ComplexF64}(undef, wave.nk)
        wr = Vector{ComplexF64}(undef, wave.nw)
        for iw in Base.OneTo(wave.nw)
            Threads.@threads for ik in Base.OneTo(wave.nk)
                nG_ik = wave.nG[ik]
                eGr_ik = view(eGr, wave.Gidx[ik])
                C = view(wave.components[ik], :, iw)
                ψr[ik] = ekr[ik] * sum(iG -> C[iG] * eGr_ik[iG], Base.OneTo(nG_ik))
            end
            wr[iw] = ekR ⋅ ψr
        end
    else
        ψr = Matrix{ComplexF64}(undef, wave.nk, 2)
        wr = Matrix{ComplexF64}(undef, 2, wave.nw)
        for iw in Base.OneTo(wave.nw)
            Threads.@threads for ik in Base.OneTo(wave.nk)
                nG_ik = wave.nG[ik]
                eGr_ik = view(eGr, wave.Gidx[ik])
                C_up = view(wave.components[ik], 1:nG_ik, iw)
                C_dn = view(wave.components[ik], (nG_ik+1):(2*nG_ik), iw)
                ψr[ik, 1] = ekr[ik] * sum(iG -> C_up[iG] * eGr_ik[iG], Base.OneTo(nG_ik))
                ψr[ik, 2] = ekr[ik] * sum(iG -> C_dn[iG] * eGr_ik[iG], Base.OneTo(nG_ik))
            end
            for ispin in 1:2
                ψr_ispin = view(ψr, :, ispin)
                wr[ispin, iw] = ekR ⋅ ψr_ispin
            end
        end
    end
    return wr
end
"""
    (wave::WannierPWKgrid)(r::AbstractVector{<:Real}, Rlist::AbstractVector{<:AbstractVector{<:Integer}})

Evaluate the Wannier functions stored in `wave` at position `r`, for the Wannier functions centered in each cell in `Rlist`.

Returns an `[nw, nR]` array if `wave.nspinor == 1`, or a `[2, nw, nR]` array if `wave.nspinor == 2`, where `nR` is the length of `Rlist`.
"""
function (wave::WannierPWKgrid{WT})(r::AbstractVector{<:Real},
    Rlist::AbstractVector{<:AbstractVector{<:Integer}}) where {WT}
    r = ReducedCoordinates(r)
    eGr = map(G->cispi(2 * (G ⋅ r)), wave.Glist)
    ekr = map(k->cispi(2 * (k ⋅ r)), wave.kgrid) ./ sqrt(wave.Ω)
    Rlist = collect(ReducedCoordinates{Int}, Rlist)
    nR = length(Rlist)
    ekR = [cispi(2 * (k ⋅ R)) / wave.nk for k in wave.kgrid, R in Rlist]
    if wave.nspinor == 1
        ψr = Vector{ComplexF64}(undef, wave.nk)
        wr = Matrix{ComplexF64}(undef, wave.nw, nR)
        for iw in Base.OneTo(wave.nw)
            Threads.@threads for ik in Base.OneTo(wave.nk)
                nG_ik = wave.nG[ik]
                eGr_ik = view(eGr, wave.Gidx[ik])
                C = view(wave.components[ik], :, iw)
                ψr[ik] = ekr[ik] * sum(iG -> C[iG] * eGr_ik[iG], Base.OneTo(nG_ik))
            end
            Threads.@threads for iR in Base.OneTo(nR)
                ekR_iR = view(ekR, :, iR)
                wr[iw, iR] = ekR_iR ⋅ ψr
            end
        end
    else
        ψr = Matrix{ComplexF64}(undef, wave.nk, 2)
        wr = Array{ComplexF64}(undef, 2, wave.nw, nR)
        for iw in Base.OneTo(wave.nw)
            Threads.@threads for ik in Base.OneTo(wave.nk)
                nG_ik = wave.nG[ik]
                eGr_ik = view(eGr, wave.Gidx[ik])
                C_up = view(wave.components[ik], 1:nG_ik, iw)
                C_dn = view(wave.components[ik], (nG_ik+1):(2*nG_ik), iw)
                ψr[ik, 1] = ekr[ik] * sum(iG -> C_up[iG] * eGr_ik[iG], Base.OneTo(nG_ik))
                ψr[ik, 2] = ekr[ik] * sum(iG -> C_dn[iG] * eGr_ik[iG], Base.OneTo(nG_ik))
            end
            Threads.@threads for iR in Base.OneTo(nR)
                ekR_iR = view(ekR, :, iR)
                for ispin in 1:2
                    ψr_ispin = view(ψr, :, ispin)
                    wr[ispin, iw, iR] = ekR_iR ⋅ ψr_ispin
                end
            end
        end
    end
    return wr
end
(wave::WannierPWKgrid)(sym::Symbol, args...; kwargs...) = wave(Val(sym), args...; kwargs...)
"""
	(wave::WannierPWKgrid{WT})(::Val{:spin_polarization}) -> Vector{real(WT)}

Return ψ_up^2 - ψ_dn^2, which is the expectation value of σz.
"""
function (wave::WannierPWKgrid{WT})(::Val{:spin_polarization}) where {WT}
    @assert wave.nspinor == 2 "Only spinor wave function need to analyse spin structure!"
    ψ_spin_polarization = zeros(real(WT), wave.nw)
    for ik in Base.OneTo(wave.nk)
        nG_ik = wave.nG[ik]
        ψ_ik = wave.components[ik]
        for iw in Base.OneTo(wave.nw)
            ψ_iw_ik_up = view(ψ_ik, 1:nG_ik, iw)
            ψ_iw_ik_dn = view(ψ_ik, (nG_ik+1):(2*nG_ik), iw)
            ψ_spin_polarization[iw] += (sum(abs2, ψ_iw_ik_up) - sum(abs2, ψ_iw_ik_dn))
        end
    end
    ψ_spin_polarization ./= wave.nk
    return ψ_spin_polarization
end
"""
	(wave::WannierPWKgrid)(::Val{:spin_index}; atol::Real = 1e-6) -> (upindex, dnindex)

If 1 - |ψ_up^2 - ψ_dn^2| < atol, then ψ obtain a spin index.
"""
function (wave::WannierPWKgrid)(::Val{:spin_index}; atol::Real=1e-6)
    @assert wave.nspinor == 2 "Only spinor wave function need to analyse spin structure!"

    ψ_spin_polarization = wave(Val(:spin_polarization))

    upindex = Vector{Int}(undef, 0)
    dnindex = Vector{Int}(undef, 0)
    spinindex = (upindex, dnindex)

    for iw in Base.OneTo(wave.nw)
        if !signbit(ψ_spin_polarization[iw])
            ispin = 1
        else
            ispin = 2
        end
        if (1 - abs(ψ_spin_polarization[iw])) < atol
            push!(spinindex[ispin], iw)
        end
    end
    return spinindex
end
"""
    translate!(wave::WannierPWKgrid, centres::Pair{<:Integer, <:AbstractVector{<:Integer}}...)

Translate the specified Wannier functions so that their centres are placed in the `[0, 0, 0]` cell.

Only specify centres `R` that are not in the `[0, 0, 0]` cell. 
This function directly modifies the `components` in `wave` so that `∑_k ψ′_{ik}` corresponds to the `[0, 0, 0]` cell.

# Example
If Wannier function 2, currently defined in the `[0, 0, 0]` cell, has its centre at `[1, 0, 0]`, call

```julia
translate!(wave, 2 => [1, 0, 0])
```
"""
function translate!(wave::WannierPWKgrid{WT},
    centres::Pair{<:Integer,<:AbstractVector{<:Integer}}...) where {WT}
    for (iw, R) in centres
        for ik in Base.OneTo(wave.nk)
            k = wave.kgrid[ik]
            enkR = cispi(-2 * (k ⋅ R))
            ψ_ik_iw = view(wave.components[ik], :, iw)
            ψ_ik_iw .*= enkR
        end
    end
    return wave
end
function Base.union(wave0::WannierPWKgrid, waves::WannierPWKgrid...)
    @info "Strict requirements: Same kgrid, Glist and Gidx, etc.!"
    waves = [wave0; waves...]
    nk = waves[1].nk
    nw = sum(wave -> wave.nw, waves)
    nspinor = waves[1].nspinor
    nG = waves[1].nG
    @assert all(wave -> wave.nk == nk, waves) "All WannierPWKgrid must have the same nk!"
    @assert all(wave -> wave.nspinor == nspinor, waves) "All WannierPWKgrid must have the same nspinor!"
    @assert all(wave -> wave.nG == nG, waves) "All WannierPWKgrid must have the same nG!"
    kgrid = waves[1].kgrid
    @assert all(wave -> wave.kgrid == kgrid, waves) "All WannierPWKgrid must have the same kgrid!"
    Gidx = waves[1].Gidx
    @assert all(wave -> wave.Gidx == Gidx, waves) "All WannierPWKgrid must have the same Gidx!"
    Glist = waves[1].Glist
    @assert all(wave -> wave.Glist == Glist, waves) "All WannierPWKgrid must have the same Glist!"
    WT = reduce(promote_type, eltype(eltype(wave.components)) for wave in waves)
    components = Vector{Matrix{WT}}(undef, nk)
    for ik in Base.OneTo(nk)
        components_ik = Matrix{WT}(undef, nG[ik], nw)
        nw_count = 0
        for wave in waves
            nw_wave = wave.nw
            components_ik[:, (nw_count+1):(nw_count+nw_wave)] .= wave.components[ik]
            nw_count += nw_wave
        end
        components[ik] = components_ik
    end
    period = waves[1].period
    @assert all(wave -> wave.period == period, waves) "All WannierPWKgrid must have the same period!"

    return WannierPWKgrid{WT}(
        nk, nw, nspinor, nG, kgrid, Gidx, Glist, components, waves[1].fft_grid, waves[1].Gmap,
        period, waves[1].lattice, waves[1].rlattice, waves[1].ldot, waves[1].rldot, waves[1].Ω, waves[1].rΩ,
    )
end
function wannier_spin_projected(wave::WannierPWKgrid{WT}; pair::Union{AbstractVector{NTuple{2,<:Integer}},Nothing}=nothing) where {WT}
    #TODO 理论上，目标子空间存在一组基全为Sz的本征态，充要条件为该子空间为Sz的不变子空间，至少，不是任意体系都可以做到。
    #TODO 这里可以应用optim方法，目标函数可以为每个wannier的(up^2-dn^2)^2，求该函数最大值，可以转化为最小值问题。
    #TODO 为了不破坏局域性，只能提供一个幺正变换矩阵，而不是每个ik一个，也就是等价于对实空间的基进行线性组合。
    #TODO 无需变换至实空间，但是需要明确wannier基属于哪个元胞，可以使用translate函数进行预操作。
    #TODO 参数为每个k点的幺正矩阵，初猜需要用户提供，也就是用户认为可以怎么变换一下，也可以直接选择恒等矩阵。
    return newwave, U
end
