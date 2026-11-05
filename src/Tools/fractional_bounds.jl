"""
    fractional_bounds(A, T_max; atol=sqrt(eps(Float64))) -> NTuple{3,Int}

Compute integer upper bounds for fractional coordinates given a Cartesian distance limit.

# Arguments
- `A`: A 3×3 matrix whose columns are basis vectors, or an `AbstractLattice` wrapping such a matrix.
- `T_max`: Maximum Cartesian distance (non-negative).
- `atol`: Absolute tolerance added before flooring to avoid missing boundary points due to floating-point round-off. Defaults to `sqrt(eps(Float64))`.

# Returns
A tuple `(N1, N2, N3)` of non-negative integers. For each fractional coordinate component `q_i`, the search range is `-N_i:N_i`.

# Mathematical background
Let `r = A * q` be the Cartesian coordinate and `G = A' * A` be the Gram matrix. The distance constraint is

    q' * G * q ≤ T_max^2.

The continuous upper bound for the i-th fractional coordinate is

    |q_i| ≤ T_max * sqrt((G^{-1})[i, i]).

The integer bound is therefore

    N_i = floor(Int, T_max * sqrt((G^{-1})[i, i]) + atol).

# Notes
The returned bounds define an axis-aligned bounding box of the ellipsoid. After enumerating integer fractional coordinates within this box, filter them with

    norm(A * q) ≤ T_max + tol

to obtain points that satisfy the distance constraint exactly (up to tolerance).
"""
function fractional_bounds(A::Union{AbstractLattice,AbstractMatrix{<:Real}}, T_max::Real; atol::Real=1e-6)
    # Unwrap lattice types to access the underlying matrix.
    A = parent(A)

    # Validate dimensions and parameters.
    size(A) == (3, 3) || throw(ArgumentError("The dimension of basis should be 3×3"))
    T_max ≥ 0 || throw(ArgumentError("T_max should be non-negative"))
    atol ≥ 0 || throw(ArgumentError("atol should be non-negative"))

    # Gram matrix G = A' * A. Convert to floating point for Cholesky.
    G = Symmetric(float.(A' * A))
    isposdef(G) || throw(ArgumentError("Basis vectors must be linearly independent"))

    # Cholesky factorization of G.
    C = cholesky(G)

    # For each axis i, solve G * x = e_i, so x = G^{-1} * e_i.
    # Then x[i] = (G^{-1})[i, i].
    N = ntuple(3) do i
        e = zeros(3)
        e[i] = 1.0
        x = C \ e
        # Guard against tiny negative values due to round-off.
        L_i = T_max * sqrt(max(x[i], zero(x[i])))
        floor(Int, L_i + atol)
    end

    return N
end
