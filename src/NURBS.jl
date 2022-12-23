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

include("bSplines/basis.jl")
include("bSplines/basisDerivatives.jl")
include("bSplines/curves.jl")
include("bSplines/surfaces.jl")

include("nurbs/basis.jl")
include("nurbs/curves.jl")
include("nurbs/surfaces.jl")

include("utils.jl")

include("utils/fileio.jl")
include("utils/plotting.jl")


end
