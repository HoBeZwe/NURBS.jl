module NURBS

# -------- used packages
using LinearAlgebra
using StaticArrays
using Requires



# -------- exportet parts
# types
#export 

# functions
export bsplineNaive
export findspan, basisfun, curvePoints
export surfacePoints, plotSurface



# -------- included files
include("bsplines.jl")
include("plotting.jl")


end
