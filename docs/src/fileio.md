
# File I/O

Two means are provided to read data from files either via the function `load` or via file specific functions.

!!! note
    To use the unifying function `load`, the package [FileIO](https://github.com/JuliaIO/FileIO.jl) is needed.



## Step files

To read [.step and .stp files](https://en.wikipedia.org/wiki/ISO_10303-21) containing B-spline surfaces and/or NURBS surfaces the `load` function and the [`readStep`](@ref readStep) function is provided.

!!! note
    So far only reading B-spline surfaces (called `B_SPLINE_SURFACE_WITH_KNOTS` in .step and .stp) and NURBS surfaces (as a `BOUNDED_SURFACE()` in .step and .stp) is supported. 
    However, reading curves should not be too dificult to implement.

```@example utils
using NURBS, FileIO
Patches = load("assets/torus.stp")

using PlotlyJS
plotPatches(Patches, plotControlPoints=false, resolution=0.25)
t = plotPatches(Patches, plotControlPoints=false, resolution=0.25) # hide
savefig(t, "torus.html"); nothing # hide
```

```@raw html
<object data="../torus.html" type="text/html"  style="width:100%;height:50vh;"> </object>
```


## Multipatch Files

To read multipatch files as defined and provided by the [nurbs Octave implementation](https://octave.sourceforge.io/nurbs/overview.html) and [BEMBEL](https://temf.github.io/bembel/) the `load` function and the [`readMultipatch`](@ref readMultipatch) function is provided.

!!! note
    Apart from these files the implementation in this package is done independendly from the Octave package. Similarities are accidentally (as it is also based on [[1]](@ref refs)). 
    Also note, that the Octave package so far provides more functionality.

```@example utils
using NURBS, FileIO
Patches = load("assets/sphere.dat")

using PlotlyJS
plotPatches(Patches, plotControlPoints=true, resolution=0.1)
t = plotPatches(Patches, plotControlPoints=true, resolution=0.1) # hide
savefig(t, "sphere.html"); nothing # hide
```

```@raw html
<object data="../sphere.html" type="text/html"  style="width:100%;height:50vh;"> </object>
```


## VTK Files (ParaView)

To display surfaces in [ParaView](https://www.paraview.org/) the function [`saveVtk`](@ref saveVtk) is provided. 
It converts the surface to an unstructured grid and exports a .vtu file which can be opened by ParaView.

```@example vtk
using NURBS, FileIO
Patches = load("assets/torus.stp")

saveVtk("filename", Patches; resolution=0.01) # filename without extension
nothing # hide
```

```@raw html
<div align="center">
<img src="../assets/paraview.png" width="750"/>
</div>
<br/>
```



To add color information the function [`vtk`](@ref NURBS.vtk) can be used together with the [WriteVKT.jl](https://juliavtk.github.io/WriteVTK.jl/stable/) package.
To this end, the color values at the four corner points of each cell have to be provided in a Vector.

```@example vtk
using NURBS, FileIO, WriteVTK
Patches = load("assets/torus.stp")

cellV, x, y, z = NURBS.vtk(Patches, 0.01)

vtk_grid("filename", x, y, z, cellV) do vtk
    vtk["dataName"] = 2 * x + y
end
nothing # hide
```

```@raw html
<div align="center">
<img src="../assets/paraviewCol.png" width="750"/>
</div>
<br/>
```

