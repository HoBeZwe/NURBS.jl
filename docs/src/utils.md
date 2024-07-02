
# Utils

For convenince the following utils are provided.


## Jacobian

To compute the [Jacobian matrix](@ref jac) and its generalized determinant a function [`Jacobian`](@ref Jacobian) is provided.

```@example utils
using NURBS # hide
Patches = readMultipatch("assets/sphere.dat")

uEvalpoints = collect(0:0.01:1.0)
vEvalpoints = collect(0:0.01:1.0)

J, dJ = Jacobian(Patches[1], uEvalpoints, vEvalpoints)
nothing # hide
```

## Greville Sites

To compute the Greville sites [[3]](@ref refs) 
```math
\gamma_i = \cfrac{u_{i+1}+ \dots + u_{i+p}}{p} \qquad i = 1, \dots, N
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