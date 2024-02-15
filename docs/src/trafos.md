
# Transformations

The following transformations can be applied to curves and surfaces.


## Scaling

To scale a curve or a surface, the functions [`scale`](@ref scale) and [`scale!`](@ref scale!) are provided.

```@example trafos
using NURBS
Patches = readMultipatch("assets/sphere.dat")

scale!(Patches, 2.3) # the sphere of radius 1 has now radius 2.3
nothing # hide
```


## Translating

To translate a curve or a surface, the functions [`translate`](@ref translate) and [`translate!`](@ref translate!) are provided.

```@example trafos
using NURBS, StaticArrays

shift = SVector(2.4,-1.0,0.5)
translate!(Patches, shift) # move the whole sphere
nothing # hide
```


## Rotating

To rotate a curve or a surface around a rotation axis by a given angle, the functions [`rotate`](@ref rotate) and [`rotate!`](@ref rotate!) are provided.

!!! tip
    The rotation axis is normalized by the rotate functions.

```@example trafos
using NURBS, StaticArrays

rotAxis = SVector(1.0, 1.0, 1.0)
α = π/3

rotate!(Patches, rotAxis, α)
nothing # hide
```


## Mirroring

To mirror/reflect a shape through a plane defined by its normal vector and an anchor point, the functions [`mirror`](@ref mirror) and [`mirror!`](@ref mirror!) are provided.

!!! tip
    The normal vector is normalized by the mirror functions.

```@example trafos
using NURBS, StaticArrays

normal = SVector(0.0, 1.0, 0.0)
anchor = SVector(0.0, 0.0, 0.0)

mirror!(Patches, normal, anchor) # mirror through the xz-plane
nothing # hide
```
