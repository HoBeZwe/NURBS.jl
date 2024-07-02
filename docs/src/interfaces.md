
# Patch Interfaces

Given a set of patches defining a surface it is often of interest to determine whether patches share common points.
To this end, functionality is provided based on the following conventions.

!!! note
    Only the cases where two patches share no points, one corner point, or one common edge are considered. Generalizations to, e.g., two common edges are possible but currently not implemented.



---
## [Per Patch Conventions](@id localNum)

On each patch the corner points and edges are assigned a (local) number as follows:
- starting at ``u=v=0``, the corners are numbered clockwise
- the edges are also numbered clockwise
- the orientations of the edges are aligned with the ``u`` and ``v`` axes

```@raw html
<div align="center">
<img src="../assets/RefSquare.svg" width="400"/>
</div>
<br/>
```


---
## Interface Conventions

When two patches share an edge, the latter is identified as interface and information is extracted about
- the IDs of the two patches
- the local edge IDs of the two patches
- the orientation of the local edges (same/opposite)
- the orientation of the edges normal to the interface (same/opposite)

When two patches share a corner point, information is exctracted about
- the IDs of the two patches
- the local corner IDs of the two patches
