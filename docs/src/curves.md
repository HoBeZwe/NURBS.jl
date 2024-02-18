
# Curves

Two type of curves are available:
- B-splines curves
- NURBS curves

Both are defined by initializing corresponding structures.

---
## Define Structures

NURBS and B-spline curves are defined by initializing a [`BsplineCurve`](@ref BsplineCurve) or a [`NURBScurve`](@ref NURBScurve) structure, respectively.

```@example curves
using NURBS # hide
# --- parameters
kVec = Float64[0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 6, 6, 6] # knot vector
kVec ./= maximum(kVec)                                # normalize it

p = 3  # degree 

# --- control points
using StaticArrays
P1 = SVector(0.0, 0.0, 0.0)
P2 = SVector(0.1, 0.25, 0.0)
P3 = SVector(0.25, 0.3, 0.0)
P4 = SVector(0.3, 0.5, 0.0)
P5 = SVector(0.4, 0.4, 0.0)
P6 = SVector(0.6, 0.3, 0.0)
P7 = SVector(0.8, 0.7, 1.0)
P8 = SVector(1.0, 0.4, 0.0)
P9 = SVector(1.1, 0.4, 0.0)

controlPoints = [P1 , P2 , P3 , P4 , P5 , P6 , P7 , P8 , P9 ]
w             = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 3.0, 1.0, 1.0] # weights for the NURBS

# --- the structures
BCurve = BsplineCurve(Bspline(p, kVec), controlPoints)
NCurve   = NURBScurve(NURB(p, kVec, w), controlPoints)
nothing # hide
```

---
## Evaluate Points of a Curve

To evaluate the curves at parametric points simply hand over the latter. 

```@example curves
evalpoints = collect(0:0.0005:1.0)

CBspline = BCurve(evalpoints)
CNurbs   = NCurve(evalpoints)
nothing # hide
```

To plot the curves the `plotCurve3D` and the `plotCurve` functions are provided, where the latter ignores any ``z``-components.

!!! note
    The [PlotlyJS.jl](https://github.com/JuliaPlots/PlotlyJS.jl) package has to be loaded in order to make the functions available.
    (It is a [weak dependency](https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions)).)

```@example curves
using PlotlyJS

plotCurve3D(CBspline, controlPoints=controlPoints)
t = plotCurve3D(CBspline, controlPoints=controlPoints) # hide
savefig(t, "cruve3D.html"); nothing # hide
```

```@raw html
<object data="../cruve3D.html" type="text/html"  style="width:100%;height:50vh;"> </object>
```

---
## Evaluate Derivatives of a Curve

To evaluate derivatives of curves at parametric points hand over the latter and the maximum derivative to be evaluated.

```@example curves
evalpoints = collect(0:0.0005:1.0)

C = NCurve(evalpoints, 2) # 0-th, 1st, and 2nd derivatives
nothing # hide
```

The `plotCurve3D` function has an optional argument `tangents` to plot vectors at the points of the curve.

```@example curves
using PlotlyJS

plotCurve3D(C[1], controlPoints=controlPoints, tangents=C[2])
t = plotCurve3D(C[1], controlPoints=controlPoints, tangents=C[2]) # hide
savefig(t, "cruve3Dder.html"); nothing # hide
```

```@raw html
<object data="../cruve3Dder.html" type="text/html"  style="width:100%;height:50vh;"> </object>
```

---
## Refining a Curve

Based on the principles of [knot insertion](@ref knotInsert), a single knot can be inserted multiple times into a curve (without changing the points the curve describes) by the [`insertKnot`](@ref insertKnot) function.

```@example curves
NCurve2 = insertKnot(NCurve, 0.75, 2) # insert knot at 0.75 twice
C2 = NCurve2(evalpoints)

plotCurve3D(C2, controlPoints=NCurve2.controlPoints)
t = plotCurve3D(C2, controlPoints=NCurve2.controlPoints) # hide
savefig(t, "cruve3DInserted.html"); nothing # hide
```

```@raw html
<object data="../cruve3DInserted.html" type="text/html"  style="width:100%;height:50vh;"> </object>
```

Several knots can be inserted by the [`refine`](@ref refine) function.

```@example curves
NCurve2 = refine(NCurve, [0.1, 0.2, 0.3, 0.8221]) # insert the array of knots
C2 = NCurve2(evalpoints)

plotCurve3D(C2, controlPoints=NCurve2.controlPoints)
t = plotCurve3D(C2, controlPoints=NCurve2.controlPoints) # hide
savefig(t, "cruve3Drefined.html"); nothing # hide
```

```@raw html
<object data="../cruve3Drefined.html" type="text/html"  style="width:100%;height:50vh;"> </object>
```

!!! note
    Refining a curve does not change the points in space described by the curve. Effectively, in the plots it can be seen that the number of control points is increased. 
    However, also underlying properties such as the differentiability are changed.


---
## Splitting a Curve

To split a curve into multiple separate curves the function [`split`](@ref NURBS.split) is provided which returns an array of curves.

```@example curves
cVec = split(NCurve, [0.2, 0.5]) # split at 0.2 and 0.5

# plot all three curves
data = PlotlyJS.GenericTrace[]
for (i, spC) in enumerate(cVec)
    push!(data, plotCurve3D(spC(evalpoints), returnTrace=true, controlPoints=spC.controlPoints)...)
end
PlotlyJS.plot(data)
t = PlotlyJS.plot(data) # hide
savefig(t, "cruve3Dsplit.html"); nothing # hide
```

```@raw html
<object data="../cruve3Dsplit.html" type="text/html"  style="width:100%;height:50vh;"> </object>
```

To equally split a curve into ``n`` curves as a second argument an integer can be provided:

```@example curves
cVec = split(NCurve, 4) # split into 4 curves
```


---
## Removing Knots from a Curve

Based on the principles of [knot removal](@ref knotRemoval), a single (interior) knot can potentially be removed multiple times from a curve (without changing the points the curve describes) by the [`removeKnot`](@ref removeKnot) function.

```@example curves
p = 3

P1 = SVector(0.0, 0.0, 1.0)
P2 = SVector(0.0, 2.0, 0.0)
P3 = SVector(1.5, 3.0, 0.0)
P4 = SVector(3.0, 3.0, 0.0)
P5 = SVector(4.5, 3.0, 0.0)
P6 = SVector(6.0, 2.0, 0.0)
P7 = SVector(6.0, 0.0, 1.0)

cP = [P1, P2, P3, P4, P5, P6, P7]

kVec = Float64[0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2]

BCurve1 = BsplineCurve(Bspline(p, kVec), cP)
C1 = BCurve1(evalpoints)

BCurve2 = removeKnot(BCurve1, 0.5, 2) # remove knot at 0.5 twice
C2 = BCurve2(evalpoints)

# --- plot both curves 
t1 = plotCurve3D(C1, controlPoints=BCurve1.controlPoints, returnTrace=true)
t2 = plotCurve3D(C2, controlPoints=BCurve2.controlPoints, returnTrace=true)

fig = make_subplots(
    rows=1, cols=2,
    specs=fill(Spec(kind="scene"), 1, 2)
)

add_trace!(fig, t1[1], row=1, col=1)
add_trace!(fig, t1[2], row=1, col=1)
add_trace!(fig, t2[1], row=1, col=2)
add_trace!(fig, t2[2], row=1, col=2)
fig
savefig(fig, "cruve3DRemoved.html"); nothing # hide
```

```@raw html
<object data="cruve3DRemoved.html" type="text/html"  style="width:100%;height:50vh;"> </object>
```

!!! note
    Removing a knot from a curve is only possible when the continuity of the curve is sufficient at the knot.
    A central part of the [`removeKnot`](@ref removeKnot) function is to verify if the knot can actually be removed.
    If not, warnings are generated, indicating the encountered limitations.