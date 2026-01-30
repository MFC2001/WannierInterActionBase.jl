
const Mat3{T} = SMatrix{3, 3, T, 9}
const Vec3{T} = SVector{3, T}

include("constant.jl")

include("./CrystallographyCore/CrystallographyCore.jl")
include("./FileFormat/FileFormat.jl")
include("./Hopping.jl")

include("./AbstractReciprocalHoppings/AbstractReciprocalHoppings.jl")

include("./SymOp.jl")
include("./BrillouinZone/BrillouinZone.jl")

include("./WilsonLoop.jl")
