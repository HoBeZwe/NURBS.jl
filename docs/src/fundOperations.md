
# Fundamental Operations

Fundamental operations to manipulate curves and surfaces are provided.

---
## Knot Insertion

To [insert a knot](@ref knotInsert) ``r`` times into a curve without changing the curve the [`insertKnot`](@ref insertKnot) function is provided.
The knot value is allowed to be already in the knot vector.

As an example consider the original curve:
```@example fundamentalOp
using NURBS # hide
# --- define curve
using StaticArrays

kVec = Float64[0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5]
kVec ./= maximum(kVec)

p = 3 # polynomial degree

P1 = SVector(0.0, 0.0, 0.0)
P2 = SVector(0.1, 0.25, 0.0)
P3 = SVector(0.25, 0.3, 0.0)
P4 = SVector(0.3, 0.5, 0.0)
P5 = SVector(0.4, 0.4, 0.0)
P6 = SVector(0.6, 0.3, 0.0)
P7 = SVector(0.8, 0.7, 1.0)
P8 = SVector(1.0, 0.4, 0.0)
P8 = SVector(1.5, 0.4, 0.0)

controlPoints = [P1, P2, P3, P4, P5, P6, P7, P8]
w             = [1.0, 0.1, 1.0, 1.0, 1.0, 1.0, 3.0, 1.0]

# --- evaluate original curve
evalpoints = collect(0:0.005:1.0)

N  = BsplineCurve(Bspline(p, kVec), controlPoints)
C1 = curvePoints(N, evalpoints)

# --- plot the curve
using PlotlyJS

plotCurve3D(C1, controlPoints=controlPoints)
t = plotCurve3D(C1, controlPoints=controlPoints) # hide
savefig(t, "curveBeforeInsert.html"); nothing # hide
```

```@raw html
<object data="../curveBeforeInsert.html" type="text/html"  style="width:120%;height:50vh;"> </object>
```



Now an already existing knot is inserted twice:
```@example fundamentalOp
uNew = 3 / 5
kVecNew, ctrlPointsNew, wNew = insertKnot(kVec, controlPoints, p, uNew, 2, w)

# --- evaluate curve with inserted point
N = BsplineCurve(Bspline(p, kVecNew), ctrlPointsNew)
C2 = curvePoints(N, evalpoints)

# --- plot the curve
plotCurve3D(C2, controlPoints=ctrlPointsNew)
t = plotCurve3D(C2, controlPoints=ctrlPointsNew) # hide
savefig(t, "curveAfterInsert.html"); nothing # hide
```

```@raw html
<object data="../curveAfterInsert.html" type="text/html"  style="width:120%;height:50vh;"> </object>
```

---
## Knot Refinement

To be done.

