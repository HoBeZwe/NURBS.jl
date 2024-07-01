
# Surfaces

Based on both bases surfaces are defined.


---
## Definitions


#### B-Spline Surfaces

A tensor product surface [[1, p. 100]](@ref refs) 
```math
\bm{s}(u,v) = \sum_{i=1}^{N_u} \sum_{j=1}^{N_v} B_{i,p}(u) B_{j,q}(v) \bm{p}_{i,j}
```
is defined by introducing two knot-vectors for the B-splines of degree ``p`` and ``q``
and a second parametric value ``v \in [0, 1]``, as well as a net of constant controlpoints ``\bm{p}_{i,j} \in \mathbb{R}^3``.

#### NURBS Surfaces

Analogously, a NURBS surface
```math
\bm{s}(u,v) = \sum_{i=1}^{N_u} \sum_{j=1}^{N_v} R_{i,j}(u,v) \bm{p}_{i,j} = \cfrac{\sum_{i=1}^{N_u} \sum_{j=1}^{N_v} B_{i,p}(u) B_{j,q}(v) w_{i,j} \bm{p}_{i,j}}{\sum_{i=1}^{N_u} \sum_{j=1}^{N_v} B_{i,p}(u) B_{j,q}(v) w_{i,j}}
```
is defined.

!!! note
    ``\bm{s}(u,v): (u, v) \mapsto (x,y,z)``, that is, the surface is a mapping from the parametric space ``(u, v) \in {[0,1]}^2`` to the physical space ``\mathbb{R}^3``.



--- 
## Derivatives


#### B-Spline Surfaces

The ``m``-th derivative in ``u`` and ``m``-th derivative in ``v`` of a B-spline surface is given as [[1, p. 111]](@ref refs)
```math
\bm{s}^{(m,n)}(u,v) = \cfrac{\partial^{m+n}}{\partial^m u \partial^n v}\bm{s}(u,v) = \sum_{i=1}^{N_u} \sum_{j=1}^{N_v} B_{i,p}^{(m)}(u) B_{j,p}^{(n)}(u) \bm{p}_{i,j} \,.
```

#### NURBS Surfaces

The ``m``-th derivative in ``u`` and ``m``-th derivative in ``v`` of a NURBS surface is given as 
```math
\bm{s}^{(m,n)}(u,v) = \cfrac{\partial^{m+n}}{\partial^m u \partial^n v}\bm{s}(u,v) = \sum_{i=1}^{N_u} \sum_{j=1}^{N_v} R_{i,j}^{(m,n)}(u,v) \bm{p}_{i,j} \,.
```
It can be computed as [[1, p. 136]](@ref refs)
```math
\begin{aligned}
\bm{s}^{(m,n)}(u,v) = \cfrac{1}{w^{(0,0)}(u,v)} \Bigg( \bm{a}^{(m,n)}(u,v) \\
&- \sum_{i=1}^m \begin{pmatrix} m \\ i \end{pmatrix} w^{(i,0)}(u,v) \bm{s}^{(m-i, n)}(u,v)\\  
&- \sum_{j=1}^n \begin{pmatrix} n \\ j \end{pmatrix} w^{(0,j)}(u,v) \bm{s}^{(m, n-j)}(u,v)\\ 
&- \sum_{i=1}^m \begin{pmatrix} m \\ i \end{pmatrix} \sum_{j=1}^n \begin{pmatrix} n \\ j \end{pmatrix}   w^{(i,j)}(u,v) \bm{s}^{(m-i, n-j)}(u,v) \Bigg) 
\end{aligned}
```
with
```math
\bm{a}^{(m,n)}(u,v) = \sum_{i=1}^{N_u} \sum_{j=1}^{N_v} B_{i,p}^{(m)}(u) B_{j,q}^{(n)}(v) w_{i,j} \bm{p}_{i,j} 
```
and
```math
w^{(m,n)}(u,v) = \sum_{i=1}^{N_u} \sum_{j=1}^{N_v} B_{i,p}^{(m)}(u) B_{j,q}^{(n)}(v) w_{i,j}  \,.
```

#### [Jacobian](@id jac)

The Jacobian matrix ``\bm{J}(u,v)`` for the surface is given as
```math
\bm{J}(u,v) = \begin{bmatrix} \cfrac{\partial x}{\partial u} & \cfrac{\partial x}{\partial v} \\[3mm] \cfrac{\partial y}{\partial u} & \cfrac{\partial y}{\partial v} \\[3mm] \cfrac{\partial z}{\partial u} & \cfrac{\partial z}{\partial v} \end{bmatrix}  = \begin{bmatrix} \bm{s}^{(1, 0)} & \bm{s}^{(0, 1)} \end{bmatrix}
```
and the magnitude of its determinant as
```math
\left|\det\left(\bm{J}(u,v)\right)\right| = \sqrt{ {\left( \cfrac{\partial y}{\partial u}\cfrac{\partial z}{\partial v} - \cfrac{\partial z}{\partial u} \cfrac{\partial y}{\partial v} \right)}^2 + {\left( \cfrac{\partial z}{\partial u}\cfrac{\partial x}{\partial v} - \cfrac{\partial x}{\partial u} \cfrac{\partial z}{\partial v} \right)}^2 + {\left( \cfrac{\partial x}{\partial u}\cfrac{\partial y}{\partial v} - \cfrac{\partial y}{\partial u} \cfrac{\partial x}{\partial v} \right)}^2 } \,.
```