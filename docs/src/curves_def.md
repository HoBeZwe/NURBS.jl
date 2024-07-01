
# Curves

Based on both bases curves are defined.

---
## Definitions

#### B-Spline Curves

Based on the B-splines ``B_{i,p}`` the ``p``-th degree curve [[1, p. 80]](@ref refs)
```math
\bm{c}(u) = \sum_{i=1}^N B_{i,p}(u) \bm{p}_i
```
is defined with the constant controlpoints ``\bm{p}_i \in \mathbb{R}^3``.


#### NURBS Curves

Analogously, the ``p``-th degree NURBS curve [[1, p. 117]](@ref refs).
```math
\bm{c}(u) = \sum_{i=1}^N R_{i,p}(u) \bm{p}_i = \cfrac{\sum_{i=1}^N B_{i,p}(u) w_i \bm{p_i}}{\sum_{i=1}^N B_{i,p}(u) w_i}
```
is defined.

!!! note
    ``\bm{c}(u): u \mapsto (x,y,z)``, that is, the curve is a mapping from the parametric space ``u \in [0,1]`` to the physical space ``\mathbb{R}^3``.


---
## Derivatives

#### B-Spline Curves

The ``k``-th derivative of a B-spline curve is given as [[1, p. 91]](@ref refs)
```math
\bm{c}^{(k)}(u) = \sum_{i=1}^N B_{i,p}^{(k)}(u) \bm{p}_i \,.
```


#### NURBS Curves

The ``k``-th derivative of a NURBS curve is given as
```math
\bm{c}^{(k)}(u) = \sum_{i=1}^N R_{i,p}^{(k)}(u) \bm{p}_i \,.
```
It can be computed as [[1, p. 125]](@ref refs)
```math
\bm{c}^{(k)}(u) = \cfrac{\bm{a}^{(k)}(u) - \sum_{i=1}^k \begin{pmatrix} k \\ i \end{pmatrix} w^{(i)}(u) \bm{c}^{(k-i)}(u)  }{w^{(0)}(u)}
```
with the auxiliary
```math
\bm{a}^{(k)}(u) = \sum_{i=1}^N B_{i,p}^{(k)}(u) w_i \bm{p}_i
```
and
```math
w^{(k)}(u) = \sum_{i=1}^N B_{i,p}^{(k)}(u) w_i \,.
```