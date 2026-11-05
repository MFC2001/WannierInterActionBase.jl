struct PlaneWaveCoulombSphere <: PlaneWaveCoulomb
    lattice::Lattice{Float64}
    Glist::Vector{ReducedCoordinates{Int}}
    R::Float64
    coulomb::Vector{Float64}
    coff::Float64
end
(pwc::PlaneWaveCoulombSphere)() = pwc.coulomb
function (pwc::PlaneWaveCoulombSphere)(::Val{:head}, qgrid::Union{AbstractVector,Tuple})
    average = pwc.R * pwc.R / 2
    return pwc.coff * average
end
"""
    PlaneWaveCoulomb(:sphere, lattice::Union{Lattice,AbstractMatrix{<:Real}}, Glist; R::Union{Real,Nothing}=nothing)

Construct a PlaneWaveCoulomb object for 1D systems with sphere truncation.
`R` is the radius of the sphere, if not provided, it will be set to half of the minimum lattice vector length.
"""
function PlaneWaveCoulomb(::Val{:sphere}, lattice::Union{Lattice,AbstractMatrix{<:Real}}, Glist;
    R::Union{Real,Nothing}=nothing)
    check_lattice_period(lattice, [false, false, false])
    if isnothing(R)
        h₁ = lattice[1, 1]
        h₂ = lattice[2, 2]
        h₃ = lattice[3, 3]
        R = min(h₁, h₂, h₃) / 2
    end

    lattice = Lattice(parent(lattice))
    rlattice = reciprocal(lattice)

    Ω = abs(det(parent(lattice)))
    coff = CoulombScale * 4π / Ω

    rlx = rlattice[1, 1]
    rly = rlattice[2, 2]
    rlz = rlattice[3, 3]
    rlx2 = rlx * rlx
    rly2 = rly * rly
    rlz2 = rlz * rlz

    Glist = ReducedCoordinates{Int}.(Glist)
    # TODO: 需要特别处理G=0
    coulomb = map(Glist) do G
        G2 = G[1] * rlx2 * G[1] + G[2] * rly2 * G[2] + G[3] * rlz2 * G[3]
        normG = sqrt(G2)
        coff * (1.0 - cos(normG * R)) / G2
    end
    return PlaneWaveCoulombSphere(lattice, Glist, R, coulomb, coff)
end
