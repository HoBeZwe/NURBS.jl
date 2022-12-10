module NURBS

# -------- used packages
using LinearAlgebra
using StaticArrays
using Requires



# -------- exportet parts
# types
export Bspline, NURB
export BsplineCurve, NURBScurve, BsplineSurface, NURBSsurface

# functions
export bSplineNaive, bSplineNaiveDerivative
export bSpline, bSplineDerivatives
export nurbsNaive

export surfacePoints, curvePoints
export readMultipatch
export generateKnotVec, numBasisFunctions



# -------- included files
include("types.jl")

include("bsplines.jl")
include("bsplineDerivatives.jl")
include("nurbs.jl")

include("utils/fileio.jl")
include("utils/plotting.jl")


end
