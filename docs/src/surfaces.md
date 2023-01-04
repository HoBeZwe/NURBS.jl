
# Surfaces

Two type of surfaces are available:
- B-splines surfaces
- NURBS surfaces

Both are defined by initializing corresponding structures.

---
## Define Structures

NURBS and B-spline surfaces are defined by initializing a [`BsplineSurface`](@ref BsplineSurface) or a [`NURBSsurface`](@ref NURBSsurface) structure, respectively.

```@example surfaces
using NURBS # hide
using StaticArrays

# --- parameters
kVec = Float64[0, 0, 0, 1, 2, 3, 4, 5, 5, 5] # knot vector
kVec ./= maximum(kVec)                       # normalize it

p = 2 # degree of the basis

# --- control points
controlPoints = [[SVector(0.0,0.0,3.0), SVector(0.0,1.0,3.0), SVector(0.0,2.0,2.0), SVector(0.0,3.0,2.0), SVector(0.0,4.0,2.0), SVector(0.0,5.0,2.0), SVector(0.0,6.0,2.0)],
                 [SVector(1.0,0.0,3.0), SVector(1.0,1.0,3.0), SVector(1.0,2.0,2.0), SVector(1.0,3.0,2.0), SVector(1.0,4.0,2.0), SVector(1.0,5.0,2.0), SVector(1.0,6.0,2.0)],
                 [SVector(2.0,0.0,2.0), SVector(2.0,1.0,2.0), SVector(2.0,2.0,1.0), SVector(2.0,3.0,1.0), SVector(2.0,4.0,1.0), SVector(2.0,5.0,1.0), SVector(2.0,6.0,1.0)],
                 [SVector(3.0,0.0,2.0), SVector(3.0,1.0,2.0), SVector(3.0,2.0,1.0), SVector(3.0,3.0,1.0), SVector(3.0,4.0,1.0), SVector(3.0,5.0,0.0), SVector(3.0,7.0,0.0)],
                 [SVector(4.0,0.0,1.0), SVector(4.0,1.0,1.0), SVector(4.0,2.0,0.0), SVector(4.0,3.0,0.0), SVector(4.0,4.0,1.0), SVector(4.0,5.0,0.0), SVector(4.0,6.0,0.0)],
                 [SVector(5.0,0.0,1.0), SVector(5.0,1.0,1.0), SVector(5.0,2.0,0.0), SVector(5.0,3.0,0.0), SVector(5.0,4.0,0.0), SVector(5.0,5.0,0.0), SVector(5.0,6.0,0.0)],
                 [SVector(6.0,0.0,1.0), SVector(6.0,1.0,1.0), SVector(6.0,2.0,0.0), SVector(6.0,3.0,0.0), SVector(6.0,4.0,0.0), SVector(6.0,5.0,0.0), SVector(6.0,6.0,0.0)]]

controlPoints = [controlPoints[i][j] for i in 1:7, j in 1:7]

# --- weights for the NURBS basis
w = ones(size(controlPoints)) 
w[5,5] = 2.0
w[7,2] = 0.8

# --- initialize structures (using the same basis in both parametric directions)
PatchB = BsplineSurface(Bspline(p, kVec), Bspline(p, kVec), controlPoints) 
PatchN = NURBSsurface(Bspline(p, kVec), Bspline(p, kVec), controlPoints, w)
nothing # hide
```

---
## Evaluate Points on the Surface

To evaluate the surfaces at parametric points the [`surfacePoints`](@ref surfacePoints) function is provided. 

```@example surfaces
uEvalpoints = collect(0:0.01:1.0)
vEvalpoints = collect(0:0.01:1.0)

SBspline = surfacePoints(PatchB, uEvalpoints, vEvalpoints)
SNurbs   = surfacePoints(PatchN, uEvalpoints, vEvalpoints)
nothing # hide
```

To plot the surfaces the `plotSurface` or the `plotPatches` functions are provided, where the latter evaluates the surface points itself.

!!! note
    The [PlotlyJS.jl](https://github.com/JuliaPlots/PlotlyJS.jl) package has to be loaded in order to make the functions available. Reason is the employed [Requires.jl](https://github.com/JuliaPackaging/Requires.jl) framework.

```@example surfaces
using PlotlyJS

plotSurface(SNurbs, controlPoints=controlPoints)
t = plotSurface(SNurbs, controlPoints=controlPoints, enforceRatio=false) # hide
savefig(t, "surface3D.html"); nothing # hide

# alternatively
plotPatches([PatchN], plotControlPoints=true)
```

```@raw html
<object data="../surface3D.html" type="text/html"  style="width:120%;height:50vh;"> </object>
```

---
## Evaluate Derivatives of the Surface

To evaluate derivatives of surfaces at parametric points the [`surfaceDerivativesPoints`](@ref surfaceDerivativesPoints) function is provided.

```@example surfaces
uEvalpoints = collect(0:0.01:1.0)
vEvalpoints = collect(0:0.01:1.0)

S = surfaceDerivativesPoints(PatchN, uEvalpoints, vEvalpoints, 2) # 0-th, 1st, and 2nd derivatives
nothing # hide
```

The `plotSurface` function has an optional argument `tangents` to plot vectors at the points of the curve.

```@example surfaces
using PlotlyJS

plotSurface(S[1, 1], tangents=S[2,1], controlPoints=PatchN.controlPoints, enforceRatio=false)
t = plotSurface(S[1, 1], tangents=S[2,1], controlPoints=PatchN.controlPoints, enforceRatio=false) # hide
savefig(t, "surface3Dder.html"); nothing # hide
```

```@raw html
<object data="../surface3Dder.html" type="text/html"  style="width:120%;height:50vh;"> </object>
```