module NURBS

# -------- used packages
using LinearAlgebra
using StaticArrays
using Requires



# -------- exportet parts
# types
export Bspline, NURB
export BsplineCurve, NURBScurve
export BsplineSurface, NURBSsurface, Surface

# functions
export evalNaive, evalNaiveDerivative
export Jacobian

export insertKnot!, insertKnot

export readMultipatch
export generateKnotVec, numBasisFunctions, spanRanges



# -------- included files
include("types.jl")

include("bSplines/basis.jl")
include("bSplines/basisDerivatives.jl")
include("bSplines/curves.jl")
include("bSplines/surfaces.jl")

include("nurbs/basis.jl")
include("nurbs/basisDerivatives.jl")
include("nurbs/curves.jl")
include("nurbs/surfaces.jl")

include("fundamentalOperations/knotInsertion.jl")

include("utils.jl")

include("utils/fileio.jl")
include("utils/plotting.jl")


end
