
# Bezier Mesh

For finite-element (FEM) applications of NURBS and B-splines it is often required to introduce a cell-structure on the NURBS surfaces that is largely independent of the surface description.
Especially, in isogeometric approaches B-splines are used to describe quantities such as currents or charges on the surface.
To this end, it is customary to introduce on each NURBS patch a virtual mesh induced by the knot vectors of the B-splines: the Bezier mesh.

!!! note
    The Bezier mesh is independent of the knot vectors of the B-splines used to describe the surface.


---
## Per Patch Conventions

The local edges and corner points numbering as defined in the [`patch conventions`](@ref localNum) is passed on to the Bezier cells.
Moreover, the numbering of the cells follows
- starting at patch one, there is a consecutive number over all patches
- on each patch the numbering starts at ``u=v=0`` and is increased first along the ``v`` axes

```@raw html
<div align="center">
<img src="../assets/BezierCells.svg" width="400"/>
</div>
<br/>
```


---
## Adjacency Information

Most importantly, adjacency information about the Bezier cells is extraced: for a given Bezier cell, what other Bezier cells:
- have a common edge at which local edge numbers
- have a common corner point at which local corner numbers