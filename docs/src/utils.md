
# Utils

For convenince the following utils are provided.

## File I/O

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
<object data="sphere.html" type="text/html"  style="width:120%;height:50vh;"> </object>
```

## Jacobian

To compute the [Jacobian matrix](@ref jac) and its generalized determinant a function [`jacobian`](@ref jacobian) is provided.

```@example utils
uEvalpoints = collect(0:0.01:1.0)
vEvalpoints = collect(0:0.01:1.0)

J, dJ = jacobian(Patches[1], uEvalpoints, vEvalpoints)
nothing # hide
```

