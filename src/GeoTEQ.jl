module GeoTEQ
using GeometryBasics

export Point, Vector

include("ma_helpers.jl")
include("Ma.jl")
include("SmoothSurface.jl")
include("IO.jl")
include("calculations.jl")

end # module GeoTEQ
