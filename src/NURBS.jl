module NURBS

# -------- used packages
using LinearAlgebra
using StaticArrays



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



# -------- extensions
function plotCurve3D end
function plotCurve end
function plotSurface end
function plotPatches end

export plotCurve3D, plotCurve, plotSurface, plotPatches



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

if !isdefined(Base, :get_extension)
    include("../ext/NURBSext.jl") # for backwards compatibility with julia versions below 1.9
end


end
