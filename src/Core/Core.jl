
export Mat3, Vec3

"""
	Mat3{T} = SMatrix{3, 3, T, 9}

A 3×3 static matrix with element type `T`.
"""
const Mat3{T} = SMatrix{3, 3, T, 9}
"""
	Vec3{T} = SVector{3, T}

A 3-component static vector with element type `T`.
"""
const Vec3{T} = SVector{3, T}

include("./CrystallographyCore/CrystallographyCore.jl")

include("constant.jl")

include("./PeriodicGridIndex.jl")
include("./broaden_peaks.jl")

include("./OpticalTransitionData.jl")


include("./SymOp.jl")
include("./BrillouinZone/BrillouinZone.jl")

include("./FileFormat/FileFormat.jl")
include("./Hopping/Hopping.jl")

include("./HR/HR.jl")
include("./PairIntHR/PairIntHR.jl")

include("./AbstractReciprocalHoppings/AbstractReciprocalHoppings.jl")
include("./HalfspaceConvexPolyhedron.jl")
include("./TripletList.jl")

include("./WannierPWKgrid.jl")

include("./WilsonLoop.jl")

include("./shell/shell.jl")
