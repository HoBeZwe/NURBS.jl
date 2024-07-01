
# [Knot Insertion](@id knotInsert)

The insertion of knots into a knot vector of a curve or a surface without changing the curve or the surface is a fundamental operation.
It is of use for [[1, p. 142]](@ref refs):
 - evaluating points and derivatives on curves and surfaces
 - subdividing curves and surfaces
 - adding control points in order to increase flexibility for interactive design control

---
## Inserting Knots into a Curve

#### Inserting a Single Knot

A given curve ``\bm{c}(u) = \sum_{i=1}^N B_{i,p}(u) \bm{p}_i`` defined on the knot vector
```math
U = \{u_1, \dots, u_M\}
```
can be represented equivalently by inserting ``\bar{u}`` into ``U`` at the ``k``th position resulting in the extended knot vector
```math
\bar{U} = \{\bar{u}_1=u_1, \dots, \bar{u}_k = u_k, \bar{u}_{k+1}=\bar{u}, \bar{u}_{k+2}=u_{k+1}, \dots, \bar{u}_{\bar{M}} = u_M \}
```
as
```math
\bm{c}(u) = \sum_{i=1}^{\bar{N}} \bar{B}_{i,p}(u) \bar{\bm{p}}_i \,,
```
where ``\bar{M} = M+1`` and ``\bar{N} = N+1`` and the ``\bar{B}_{i,p}`` are the ``p``th degree basis functions defined on ``\bar{U}``.

!!! note
    Hence, there are two tasks to fulfill:
    - determining the position index ``k`` and the resulting ``\bar{U}`` for a given ``\bar{u}``
    - determining the ``\bar{\bm{p}}_i``

While the first one is a standard task, the second one can be computed as [[1, p. 143]](@ref refs)
```math
\bar{\bm{p}}_i = \alpha_i \bm{p}_i + (1 - \alpha_i) \bm{p}_{i-1}
```
with the coefficients
```math
\alpha_i = \begin{cases} 1 & i \leq k - p \\ \cfrac{\bar{u} - u_i }{ u_{i+p} - u_i} & k - p + 1 \leq i \leq k \\ 0 & i \geq k + 1 \end{cases} \,.
```
In consequence, ``p`` new control points have to be computed replacing ``p-1`` existing ones.


#### Inserting Mulitple Knots

If the point ``\bar{u}`` to be inserted, already exists ``s`` times in ``U`` and shall be inserted ``r`` additional times, the above formula generalizes to [[1, p. 149]](@ref refs)
```math
\bar{\bm{p}}_{i,r} = \alpha_{i,r} \bar{\bm{p}}_{i, r-1} + (1 - \alpha_{i,r}) \bar{\bm{p}}_{i-1, r-1}
```
with ``\bar{\bm{p}}_{i,0} = \bm{p}_i`` and
```math
\alpha_{i,r} = \begin{cases} 1 & i \leq k - p + r - 1 \\ \cfrac{\bar{u} - u_i }{ u_{i+p-r+1} - u_i} & k - p + r \leq i \leq k - s \\ 0 & i \geq k - s + 1 \end{cases} \,.
```

!!! note
    - ``p-s+r-1`` new control points are computed.
    - ``p-s-1`` control points are replaced starting at index ``k - p + 1``.


---
## Inserting Knots into a Surface

Inserting a knot into a surface (either along ``u`` or along ``v``) basically corresponds to applying the knot insertion process from a curve to the grid of control points of the surface.


---
## Knot Refinement

Knot refinement is the insertion of several different knots into a surface or a curve.


---
## Splitting Curves and Surfaces

By inserting a knot into a surface such that the knot multiplicity is equal to the degree of the curve, the curve can be separated into two curves.
The two curves describe the same points in space as the original curve.

Analogously, a surface can be split into multiple surfaces via knot insertion.
