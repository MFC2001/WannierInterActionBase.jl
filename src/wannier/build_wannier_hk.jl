"""
	build_wannier_hk(bands, chkdata::wannier90_chk) -> Vector{Matrix}
	build_wannier_hk(bands, u_matrix_total; bandindex=nothing) -> Vector{Matrix}

Build the Wannier-basis Hamiltonian `H(k) = U(k)' * diag(bands[k]) * U(k)` for every
k-point and return the vector of `nw × nw` Hermitian matrices. With a
`wannier90_chk`, the bands in the window are selected automatically; otherwise
`bandindex[ik]` may restrict the bands used at each k-point.
"""
function build_wannier_hk(bands::AbstractMatrix{<:Real}, chkdata::wannier90_chk)
    return build_wannier_hk(bands, chkdata.u_matrix_total; bandindex=chkdata.inwindow_bands)
end
function build_wannier_hk(bands::AbstractMatrix{ET}, u_matrix_total::AbstractVector{<:AbstractMatrix{T}};
    bandindex::Union{AbstractVector{<:AbstractVector{<:Integer}},Nothing}=nothing) where {ET<:Real,T<:Number}

    (_, nk) = size(bands)
    @assert length(u_matrix_total) == nk "Mismatched bands and u_matrix."

    bands_vec = Vector{Vector{ET}}(undef, nk)
    if isnothing(bandindex)
        for ik in 1:nk
            bands_vec[ik] = bands[:, ik]
        end
    else
        for ik in 1:nk
            bands_vec[ik] = bands[bandindex[ik], ik]
        end
    end

    return build_wannier_hk(bands_vec, u_matrix_total)
end
function build_wannier_hk(bands::AbstractVector{<:AbstractVector{<:Real}},
    u_matrix_total::AbstractVector{<:AbstractMatrix{T}}) where {T<:Number}

    nk = length(bands)
    @assert length(u_matrix_total) == nk "Mismatched bands and u_matrix."
    nw = size(u_matrix_total[1], 2)
    HT = T <: Complex ? ComplexF64 : Float64
    Hk = Vector{Matrix{HT}}(undef, nk)

    for ik in 1:nk
        bands_ik = bands[ik]
        u_matrix_ik = u_matrix_total[ik]
        Hk_ik = Matrix{HT}(undef, nw, nw)
        for j in 1:nw, i in 1:j
            hk_ij = sum(ib -> bands_ik[ib] * conj(u_matrix_ik[ib, i]) * u_matrix_ik[ib, j], eachindex(bands_ik))
            Hk_ik[i, j] = hk_ij
            Hk_ik[j, i] = conj(hk_ij)
        end
        Hk[ik] = Hk_ik
    end

    return Hk
end
