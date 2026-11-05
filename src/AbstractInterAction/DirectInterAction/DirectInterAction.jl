
include("./LongRangePart/LongRangePart.jl")

export DirectInterAction
"""
	DirectInterAction <: AbstractInterAction

Direct Coulomb interaction defined in reciprocal space.

\$V_{ij}(\\mathbf{k})=\\sum_{\\mathbf{R}} V_{i\\mathbf{0}j\\mathbf{R}} e^{i\\mathbf{k}\\cdot\\mathbf{R}}\$
"""
abstract type DirectInterAction <: AbstractInterAction end

include("./DirectLR.jl")
include("./DirectSR.jl")
include("./DirectSRLR.jl")

numorb(U::DirectInterAction) = U.norb
(v::DirectInterAction)(sym::Symbol, args...; kwargs...) = v(Val(sym), args...; kwargs...)
# nearΓ代表计算时特别处理发散项，根据输入参数是否包含nkorkgrid，其行为分别为舍去发散项(仅舍弃发散部分，例如三维的1/k^2)和将发散项替换为头项积分值。
# isΓ在nearΓ=true的情况下仅用作强调，内部会在计算上做一些简化。
# 可以仅指定isΓ，此时nearΓ = isΓ。
# isΓ = true, nearΓ = false的情况下，等同于将k点视为正常点，这时如果k趋于Γ点，返回值不能保证一定是正常值。
# 提供方法v(:G0scale, k)，以便于获取发散项中去除发散部分后的系数，也就是k²*infinite_term(三维)或者k*infinite_term(三维)的值。
# k的周期为[-0.5, 0.5]区间，超出这个区间的k会被自动转换到这个区间内。
function (v::DirectInterAction)(k::AbstractVector{<:Real}; isΓ::Bool=false, nearΓ::Bool=isΓ)
    A = Matrix{ComplexF64}(undef, v.norb, v.norb)
    return v(A, ReducedCoordinates(k); isΓ, nearΓ)
end
function (v::DirectInterAction)(k::AbstractVector{<:Real},
    nkorkgrid::Union{Integer,AbstractVector{<:Integer},NTuple{3,<:Integer}}; isΓ::Bool=false, nearΓ::Bool=isΓ)
    A = Matrix{ComplexF64}(undef, v.norb, v.norb)
    return v(A, ReducedCoordinates(k), nkorkgrid; isΓ, nearΓ)
end

"""
	DirectInterAction(sym::Symbol, U, lattice, orblocat::AbstractVector; kwargs...) -> DirectInterAction

Construct a `DirectSRLR` instance with given long-range ineraction `sym`, interaction `U` in a finite range, lattice, orbital locations.
This instance contains long-range interaction and short-range interaction.

	DirectInterAction(sym::Symbol, lattice, orblocat::AbstractVector; kwargs...) -> DirectInterAction

Construct a `DirectLR` instance with given long-range ineraction `sym`, lattice, orbital locations.
This instance only contains long-range interaction.

	DirectInterAction(U) -> DirectInterAction

Construct a `DirectSR` instance with given interaction `U` in a finite range.
This instance only contains the interaction in `U`, and acctually it's just a simple packaging of U.

- `U`: the direct Coulomb potential between the wannier bases within a limited distance, can be a `HR` or `AbstractReciprocalHoppings`.
kwargs:
- `period`: the periodicity;
- `ϵ::Real`: dielectric constant;
- `αrcut::Real`: α × rcut, default is 4.5;
- `δ:Real`: the cutoff value of reciprocal potential.

Make sure the input U has converged to 1/r.
"""
function DirectInterAction(sym::Symbol, U::AbstractHR,
    lattice::Union{Lattice,AbstractMatrix{<:Real}}, orblocat::AbstractVector{V}; kwargs...) where {V}

    # Note: AbstractHR contains ndegen, which only can be applied to U_short, 
    # so it should be transform to ReciprocalHoppings after minus ϕ!
    numorb(U) == length(orblocat) || error("Wrong number of orbital locations.")
    orblocat = collect(orblocat)
    lattice = Lattice(parent(lattice))
    if V <: CartesianCoordinates
        orblocat = map(orblocat) do r
            lattice \ r
        end
    else
        orblocat = map(orblocat) do r
            ReducedCoordinates(r)
        end
    end
    return DirectSRLR(Val(sym), U, lattice, orblocat; kwargs...)
end
function DirectInterAction(sym::Symbol,
    lattice::Union{Lattice,AbstractMatrix{<:Real}}, orblocat::AbstractVector{V}; kwargs...) where {V}
    orblocat = collect(orblocat)
    norb = length(orblocat)
    lattice = Lattice(parent(lattice))
    if V <: CartesianCoordinates
        orblocat = map(orblocat) do r
            lattice \ r
        end
    else
        orblocat = map(orblocat) do r
            ReducedCoordinates(r)
        end
    end
    V_LR = LongRangePart(Val(sym), lattice, orblocat; kwargs...)
    return DirectLR(norb, V_LR)
end
function DirectInterAction(U::AbstractHR)
    U = ReciprocalHoppings(U)
    return DirectSR(numorb(U), U)
end
