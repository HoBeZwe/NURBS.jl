module NURBS

# -------- used packages
using LinearAlgebra
using StaticArrays
using Requires



# -------- exportet parts
# types
#export 

# functions
export bSplineNaive, bSplineNaiveDerivative
export findSpan, basisFun, derBasisFun, curvePoints
export surfacePoints, plotSurface



# -------- included files
include("bsplines.jl")
include("bsplineDerivatives.jl")
include("nurbs.jl")
include("plotting.jl")


end
