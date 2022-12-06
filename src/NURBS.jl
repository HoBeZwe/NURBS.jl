module NURBS

# -------- used packages
using LinearAlgebra
using StaticArrays
using Requires



# -------- exportet parts
# types
export Bspline, NURB
export Curve, BsplineSurface, NURBSsurface

# functions
export bSplineNaive, bSplineNaiveDerivative
export findSpan, basisFun, derBasisFun, curvePoints
export surfacePoints
export readMultipatch



# -------- included files
include("types.jl")

include("bsplines.jl")
include("bsplineDerivatives.jl")
include("nurbs.jl")

include("utils/fileio.jl")
include("utils/plotting.jl")


end
