module NURBS

# -------- used packages
using LinearAlgebra
using Statistics
using StaticArrays
using Suppressor



# -------- exportet parts
# types
export Bspline, NURB, CurrySchoenberg
export BsplineCurve, NURBScurve
export BsplineSurface, NURBSsurface, Surface

# functions
export evalNaive, evalNaiveDerivative
export Jacobian, JacobiDet

export insertKnot!, insertKnot
export refine # NOTE: split is an extension of Base.split, i.e., it is publicly available without export
export removeKnot!, removeKnot

export readMultipatch, readStep
export generateKnotVec, numBasisFunctions, spanRanges
export greville, anchors
export scale, scale!
export translate, translate!
export rotate, rotate!
export mirror, mirror!


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
include("fundamentalOperations/splitting.jl")
include("fundamentalOperations/knotRemoval.jl")

include("utils.jl")
include("fileIO/multipatch.jl")
include("fileIO/step.jl")

if !isdefined(Base, :get_extension)
    include("../ext/NURBSext.jl") # for backwards compatibility with julia versions below 1.9
end


end
