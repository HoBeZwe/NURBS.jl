
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
BSurface = BsplineSurface(Bspline(p, kVec), Bspline(p, kVec), controlPoints) 
NSurface = NURBSsurface(Bspline(p, kVec), Bspline(p, kVec), controlPoints, w)
nothing # hide
```

---
## Evaluate Points on a Surface

To evaluate the surfaces at parametric points hand over the latter. 

```@example surfaces
uEvalpoints = collect(0:0.01:1.0)
vEvalpoints = collect(0:0.01:1.0)

SBspline = BSurface(uEvalpoints, vEvalpoints)
SNurbs   = NSurface(uEvalpoints, vEvalpoints)
nothing # hide
```

To plot the surfaces the `plotSurface` or the `plotPatches` functions are provided, where the latter evaluates the surface points itself.

!!! note
    The [PlotlyJS.jl](https://github.com/JuliaPlots/PlotlyJS.jl) package has to be loaded in order to make the functions available.
    (It is a [weak dependency](https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions)).)

```@example surfaces
using PlotlyJS

plotSurface(SNurbs, controlPoints=controlPoints)
t = plotSurface(SNurbs, controlPoints=controlPoints, enforceRatio=false) # hide
savefig(t, "surface3D.html"); nothing # hide

# alternatively
#plotPatches([NSurface], plotControlPoints=true)
```

```@raw html
<object data="../surface3D.html" type="text/html"  style="width:100%;height:50vh;"> </object>
```

---
## Evaluate Derivatives of a Surface

To evaluate derivatives of surfaces at parametric points hand over the latter and the maximum derivative to be evaluated.

```@example surfaces
uEvalpoints = collect(0:0.01:1.0)
vEvalpoints = collect(0:0.01:1.0)

S = NSurface(uEvalpoints, vEvalpoints, 2) # 0-th, 1st, and 2nd derivatives
nothing # hide
```

In case points (e.g., single points) shall be evaluated many times on demand, memory can be preallocated and reused in subsequent calls:

```@example surfaces
pM = NURBS.preAllocNURBSsurface(p, p, uEvalpoints, vEvalpoints, 2)

S = NSurface(uEvalpoints, vEvalpoints, 2, pM)
nothing # hide
```


The `plotSurface` function has an optional argument `tangents` to plot vectors at the points of the curve.

```@example surfaces
using PlotlyJS

plotSurface(S[1, 1], tangents=S[2,1], controlPoints=NSurface.controlPoints, enforceRatio=false)
t = plotSurface(S[1, 1], tangents=S[2,1], controlPoints=NSurface.controlPoints, enforceRatio=false) # hide
savefig(t, "surface3Dder.html"); nothing # hide
```

```@raw html
<object data="../surface3Dder.html" type="text/html"  style="width:100%;height:50vh;"> </object>
```

---
## Refining a Surface

Based on the principles of [knot insertion](@ref knotInsert), several knots can be inserted into a surface (without changing the points the surface describes) by the [`refine`](@ref refine) function.

```@example surfaces
BSurfaceRef = refine(BSurface; U=[0.5, 0.5, 0.72], V=[0.2, 0.42])

plotPatches([BSurfaceRef], enforceRatio=false, plotControlPoints=true)
t = plotPatches([BSurfaceRef], enforceRatio=false, plotControlPoints=true) # hide
savefig(t, "surfaceRefined.html"); nothing # hide
```

```@raw html
<object data="../surfaceRefined.html" type="text/html"  style="width:100%;height:50vh;"> </object>
```

!!! note
    Refining a surface does not change the points in space described by the surface. Effectively, in the plots it can be seen that the number of control points is increased. 
    However, also underlying properties such as the differentiability are changed.


---
## Splitting a Surface

To split a surface into multiple separate surfaces the function [`split`](@ref NURBS.split) is provided which returns an array of surfaces.

```@example surfaces
sVec = split(BSurface, U=[0.5,0.75], V=[0.6])

plotPatches(sVec, enforceRatio=false, plotControlPoints=false)
t = plotPatches(sVec, enforceRatio=false, plotControlPoints=false) # hide
savefig(t, "surfaceSplit.html"); nothing # hide
```

```@raw html
<object data="../surfaceSplit.html" type="text/html"  style="width:100%;height:50vh;"> </object>
```

To equally split a surface into ``n`` times ``m`` surfaces as second and third argument an integer can be provided:

```@example surfaces
sVec = split(BSurface, 3, 2) # split into 3 x 2 surfaces
```


---
## Removing Knots from a Surface

Based on the principles of [knot removal](@ref knotRemoval), an interior knot can potentially be removed multiple times from a surface (without changing the points the surface describes) by the [`removeKnotU`](@ref removeKnotU) and the [`removeKnotV`](@ref removeKnotV) functions.

```@example surfaces
# --- surface with removable knots
p = 3
kVec = Float64[0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2]

P1 = SVector(0.0, 0.0, 0.0)
P2 = SVector(0.0, 2.0, 0.0)
P3 = SVector(1.5, 3.0, 0.0)
P4 = SVector(3.0, 3.0, 0.0)
P5 = SVector(4.5, 3.0, 0.0)
P6 = SVector(6.0, 2.0, 0.0)
P7 = SVector(6.0, 0.0, 0.0)

cP = [P1 , P2 , P3 , P4 , P5 , P6 , P7 ]

P1 = SVector(0.0, 0.0, 2.0)
P2 = SVector(0.0, 2.0, 2.0)
P3 = SVector(1.5, 3.0, 2.0)
P4 = SVector(3.0, 3.0, 2.0)
P5 = SVector(4.5, 3.0, 2.0)
P6 = SVector(6.0, 2.0, 2.0)
P7 = SVector(6.0, 0.0, 2.0)

cP2 = [P1 , P2 , P3 , P4 , P5 , P6 , P7 ]

controlPoints = [[cP, cP2][j][i] for i in 1:7, j in 1:2] 

BS1 = BsplineSurface(Bspline(p, kVec), Bspline(1, [0.0,0.0,1.0,1.0]), controlPoints)


# --- remove knot
BS2 = removeKnotU(BS1, 0.5, 1) # remove knot once


# --- plot before and after
t1 = plotPatches([BS1], enforceRatio=true, plotControlPoints=true)
t2 = plotPatches([BS2], enforceRatio=true, plotControlPoints=true)

uEval = vEval = collect(0:0.005:1.0)

S1 = BS1(uEval, vEval)
S2 = BS2(uEval, vEval)

x, t1 = plotSurface(S1, controlPoints=BS1.controlPoints, returnTrace=true)
x, t2 = plotSurface(S2, controlPoints=BS2.controlPoints, returnTrace=true)

fig = make_subplots(
    rows=1, cols=2,
    specs=fill(Spec(kind="scene"), 1, 2)
)

for i in eachindex(t1)
    add_trace!(fig, t1[i], row=1, col=1)
end
for i in eachindex(t2)
    add_trace!(fig, t2[i], row=1, col=2)
end

fig
savefig(fig, "surfaceRemoved.html"); nothing # hide
```

```@raw html
<object data="../surfaceRemoved.html" type="text/html"  style="width:100%;height:50vh;"> </object>
```

!!! note
    Removing a knot from a surface is only possible when the continuity of the surface is sufficient at the knot.
    A central part of the removeKnot functions is to verify if the knot can actually be removed.
    If not, warnings are generated, indicating the encountered limitations.