
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
NBspline = BsplineCurve(Bspline(p, kVec), controlPoints)
NNurbs   = NURBScurve(NURB(p, kVec, w), controlPoints)
nothing # hide
```

---
## Evaluate Points on the Curve

To evaluate the curves at parametric points simply hand over the latter. 

```@example curves
evalpoints = collect(0:0.0005:1.0)

CBspline = NBspline(evalpoints)
CNurbs   = NNurbs(evalpoints)
nothing # hide
```

To plot the curves the `plotCurve3D` and the `plotCurve` functions are provided, where the latter ignores any 'z'-components.

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
<object data="../cruve3D.html" type="text/html"  style="width:120%;height:50vh;"> </object>
```

---
## Evaluate Derivatives of the Curve

To evaluate derivatives of curves at parametric points hand over the latter and the maximum derivative to be evaluated.

```@example curves
evalpoints = collect(0:0.0005:1.0)

C = NNurbs(evalpoints, 2) # 0-th, 1st, and 2nd derivatives
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
<object data="../cruve3Dder.html" type="text/html"  style="width:120%;height:50vh;"> </object>
```