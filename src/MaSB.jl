module MaSB

using GeometryBasics

export Point, Vector


include("ma_helpers.jl")
include("Ma.jl")
include("SmoothSurface.jl")

end # module MaSB
