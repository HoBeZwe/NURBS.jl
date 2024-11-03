
# Connectivity

The connectivity of a set of surface patches can be identified following the [`patch conventions`](@ref localNum) and the [`Bezier mesh conventions`](@ref localNum).

---
## Patch Interfaces

To identify common edges and common corner points beween patches of a surface, the [`identifyInterfaces`](@ref identifyInterfaces) function is provided

```@example meshes
using NURBS
Patches = readMultipatch("assets/sphere.dat")

interfaces, commonVtxs = identifyInterfaces(Patches, tol=1e-3)
nothing # hide
```

!!! note 
    To determine whether two patches share a common edge, in the current implementation only the corner points are compared. The maximim allowed distance between two points to be identified as one can be specified by `tol`.

The returns are vectors of [`Interface`](@ref Interface) structures, which are a [winged edge](https://en.wikipedia.org/wiki/Winged_edge) like structure:
```@example
using StaticArrays # hide
struct PatchInterface{I}
    patchID::I
    localEdge::I
end

struct Interface{I}
    patch1::PatchInterface{I}
    patch2::PatchInterface{I}
    orientation::I              # edges normal to interface have same orientation
    reverse::Bool               # adjacent edges have same orientation
end
```

To get patch-wise information the [`getPatchInterfaces`](@ref getPatchInterfaces) function is provided
```@example meshes
ifPatchwise = NURBS.getPatchInterfaces(Patches, interfaces, commonVtxs)
nothing # hide
```



---
## Bezier Mesh

To extract information about the adjacency of a virtual Bezier mesh, the 
```@example meshes
bezierAdj = NURBS.bezierAdjacency(interfaces, commonVtxs, 3, 3, length(Patches))
nothing # hide
```

!!! tip
    To determine all connectivity information, the [`Connectivity`](@ref NURBS.Connectivity) function can be used, which returns a struct containing the interface vector, the patchwise vectors and the Bezier adjacency:

    ```@example meshes
    using NURBS

    Patches = readMultipatch("assets/sphere.dat")

    cty = Connectivity(Patches, 3, 3, tol=1e-3)

    cty.interfaces
    cty.patchInterfaces
    cty.bezierAdjacency
    nothing # hide
    ```