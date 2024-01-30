
# Utils

For convenince the following utils are provided.

## File I/O

Two means are provided to read data from files.

### Step files

To read [.step and .stp files](https://en.wikipedia.org/wiki/ISO_10303-21) containing B-spline surfaces and/or NURBS surfaces a function [`readStep`](@ref readStep) is provided.

!!! note
    So far only reading B-spline surfaces (called `B_SPLINE_SURFACE_WITH_KNOTS` in .step and .stp) and NURBS surfaces (as a `BOUNDED_SURFACE()` in .step and .stp) is supported. 
    Howver, reading curves should not be too dificult to implement.

```@example utils
using NURBS # hide
Patches = readStep("assets/torus.stp")

using PlotlyJS
plotPatches(Patches, plotControlPoints=false, resolution=0.25)
t = plotPatches(Patches, plotControlPoints=false, resolution=0.25) # hide
savefig(t, "torus.html"); nothing # hide
```

```@raw html
<object data="../torus.html" type="text/html"  style="width:100%;height:50vh;"> </object>
```


### Multipatch Files

To read multipatch files as defined and provided by the [nurbs Octave implementation](https://octave.sourceforge.io/nurbs/overview.html) and [BEMBEL](https://temf.github.io/bembel/) a function [`readMultipatch`](@ref readMultipatch) is provided.

!!! note
    Apart from these files the implementation in this package is done independendly from the Octave package. Similarities are accidentally (as it is also based on [[1]](@ref refs)).
    Also note, that the Octave package so far provides more functionality.

```@example utils
using NURBS # hide
Patches = readMultipatch("assets/sphere.dat")

using PlotlyJS
plotPatches(Patches, plotControlPoints=true, resolution=0.1)
t = plotPatches(Patches, plotControlPoints=true, resolution=0.1) # hide
savefig(t, "sphere.html"); nothing # hide
```

```@raw html
<object data="../sphere.html" type="text/html"  style="width:100%;height:50vh;"> </object>
```

## Jacobian

To compute the [Jacobian matrix](@ref jac) and its generalized determinant a function [`Jacobian`](@ref Jacobian) is provided.

```@example utils
uEvalpoints = collect(0:0.01:1.0)
vEvalpoints = collect(0:0.01:1.0)

J, dJ = Jacobian(Patches[1], uEvalpoints, vEvalpoints)
nothing # hide
```

## Greville Sites

To compute the Greville sites [[3]](@ref refs) 
```math
\gamma_i = \cfrac{u_{i+1}+ \dots + u_{i+p}}{p} \qquad i = 1, \dots, B
```
corresponding to a given knot vector with entries ``u_i`` and a polynomial degree ``p`` the function [`greville`](@ref greville) is provided.

```@example utils
p = 2

kVec = generateKnotVec(5, p)
Bspl = Bspline(p, kVec)

gs = greville(Bspl)
```

## Anchor Sites

To compute the anchors [[4]](@ref refs) corresponding to a given knot vector and a polynomial degree ``p`` the function [`anchors`](@ref anchors) is provided.

```@example utils
ac = anchors(Bspl)
```