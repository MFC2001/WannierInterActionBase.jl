"""
Some data about lammps.
"""
function headerdict()::Dict{String,Any}
    return Dict(
        "outorder" => ["atoms", "atom types",
            "xlo xhi", "ylo yhi", "zlo zhi", "xy xz yz"],
        "defaut" => 0,
        "atoms" => 0,
        "atom types" => 0,
        "xlo xhi" => [-0.5, 0.5],
        "ylo yhi" => [-0.5, 0.5],
        "zlo zhi" => [-0.5, 0.5],
        "xy xz yz" => [0, 0, 0]
    )
end
function bodydict()::Dict{String,Any}
    return Dict(
        "outorder" => ["Masses", "Atoms", "Velocities"],
        "Masses" => :(Masses),
        "Atoms" => :(Atoms),
        "Velocities" => :(Velocities)
    )
end
function atom_styledict()::Dict{String,Any}
    #Can replace empty arraies with function names if needed.
    return Dict(
        "atomic" => [],
    )
end
struct AtomsData
    index::Vector{Int}
    type::Vector{Int}
    location::Matrix{Float64}
end
struct MassesData
    type::Vector{Int}
    value::Vector{Float64}
end