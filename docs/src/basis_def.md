
# Bases

The two considered B-spline and NURBS basis functions are based on a knot vector.

---
## [Knot Vectors](@id kvec)

Given a polynomial degree ``p \in \{1, 2, 3, \dots\}``, a knot vector of length ``M`` is defined as 

```math
U = \{u_1, \dots, u_M\}\,,
```

where ``u_i \leq u_{i+1}``, i.e., the entries are sorted in ascending order, ``u_1 = 0``, ``u_M = 1``, and the first and last entry are each repeated ``p + 1`` times.


!!! note
    This specific definition of the first and last entries defines an **open** knot vector; the only type of knot vector considered in this package.

!!! note
    Having ``u_i \in [0, 1]`` is assumed (and partially enforced) everywhere. 


---
## Definitions 

### [B-Splines](@id bspl)

The B-spline basis functions of degree ``p`` are defined as ``N_{i,p}`` on the knot vector ``U`` recursively, starting with the piecewise constant (``p=0``) [[1, p. 50]](@ref refs)
```math
N_{i,0}(u) = \begin{cases} 1 & u_i \leq u < u_{i+1} \\ 0 & \text{otherwise}\end{cases}
```
as well as for ``p>0``
```math
N_{i,p}(u) = \cfrac{u - u_i}{u_{i+p} - u_i} \, N_{i, p-1}(u) + \cfrac{u_{i+p+1} - u}{u_{i+p+1} - u_{i+1}} \, N_{i+1, p-1}(u) \,.
```
Whenever one of the contained quotients exhibits a division by 0, the quotient is defined to be 0 itself.

!!! note
    The number of basis functions ``B`` is related to the length of the knot vector ``M`` and the polynomial degree ``p`` as
    ```math
    B = M - p - 1 \,.
    ```


### [NURBS](@id nurbs)

Introducing ``B`` weights ``w_i \in \mathbb{R}_+`` the rational B-spline basis functions [[1, p. 118]](@ref refs)
```math
R_{i, p}(u) = \cfrac{N_{i, p}(u) w_i}{\sum_{j=1}^B N_{j, p}(u) w_i}
```
are defined based on the B-spline basis ``N_{i, p}`` on the knot vector ``U``.


---
## [Derivatives](@id derB)

### B-Splines

The ``k``-th derivative of the B-splines ``N_{i, p}`` can be computed as [[1, p. 61]](@ref refs)
```math
N_{i,p}^{(k)}(u) = p \left( \cfrac{N_{i,p-1}^{(k-1)}}{u_{i+p} - u_i} - \cfrac{N_{i+1,p-1}^{(k-1)}}{u_{i+p+1} - u_{i+1}} \right) \,.
```

### NURBS

The ``k``-th derivative of the NURBS basis functions ``R_{i, p}`` can be computed as [[2]](@ref refs)
```math
R_{i,p}^{(k)}(u) = \cfrac{w_i \, N_{i,p}^{(k)}(u) - \sum_{j=1}^k \begin{pmatrix} k \\j\end{pmatrix} W^{(j)}(u) R_{i,p}^{k-j}(u)}{W^{(0)}(u)}
```
with
```math
W^{(k)}(u) = \sum_{i=1}^B N_{i, p}^{(k)}(u) w_i \,.
```