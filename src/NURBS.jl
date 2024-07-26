module NURBS

# -------- used packages
using LinearAlgebra
using Statistics
using StaticArrays
using Suppressor
using FileIO
using UUIDs



# -------- exportet parts
# types
export Bspline, NURB, CurrySchoenberg
export BsplineCurve, NURBScurve
export BsplineSurface, NURBSsurface, Surface
export degree

# functions
export evalNaive, evalNaiveDerivative
export Jacobian, JacobiDet

export insertKnot!, insertKnot
export coarsen, refine # NOTE: split is an extension of Base.split, i.e., it is publicly available without export
export removeKnot!, removeKnot, removeKnotU, removeKnotV

export readMultipatch, readStep
export generateKnotVec, numBasisFunctions, spanRanges
export greville, anchors
export scale, scale!
export translate, translate!
export rotate, rotate!
export mirror, mirror!

export PatchInterface, Interface, InterfacePatchwise, InterfaceData, Connectivity
export identifyInterfaces, getPatchInterfaces, patchID, localEdge, patchIDs, localEdges

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
include("fileIO/fileIO.jl")

function __init__() # locally initialize file formats for FileIO load and save functions
    add_format(format"DAT", io -> endswith(io.name, ".dat>"), [".dat"], [:NURBS => UUID("dde13934-061e-461b-aa91-2c0fad390a0d")])
    add_format(
        format"STEP",
        "ISO-10303-21",
        [".stp", ".step", ".stpnc", ".p21", ".210"],
        [:NURBS => UUID("dde13934-061e-461b-aa91-2c0fad390a0d")],
    )
end
# NOTE: think about registering reading of step files officially to the FileIO package


include("connectivity/interfaces.jl")
include("connectivity/bezierMesh.jl")

if !isdefined(Base, :get_extension)
    include("../ext/NURBSext.jl") # for backwards compatibility with julia versions below 1.9
end


end
