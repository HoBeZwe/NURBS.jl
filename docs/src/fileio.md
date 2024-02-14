
# File I/O

Two means are provided to read data from files.

## Step files

To read [.step and .stp files](https://en.wikipedia.org/wiki/ISO_10303-21) containing B-spline surfaces and/or NURBS surfaces a function [`readStep`](@ref readStep) is provided.

!!! note
    So far only reading B-spline surfaces (called `B_SPLINE_SURFACE_WITH_KNOTS` in .step and .stp) and NURBS surfaces (as a `BOUNDED_SURFACE()` in .step and .stp) is supported. 
    However, reading curves should not be too dificult to implement.

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


## Multipatch Files

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