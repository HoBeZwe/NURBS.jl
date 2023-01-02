
# NURBS.jl

This package provides B-spline and NURBS (non-uniform rational B-spline) basis functions, their derivatives, as well as curves and surfaces based on both considered basis functions.


---
## Overview

The following aspects are implemented (✓) and planned (⌛):

##### B-spline & NURBS
- ✓ Basis & derivatives
- ✓ Curves & derivatives
- ✓ Surfaces & derivatives

##### Fundamental operations
- ✓ File I/O (basic)
- ⌛ Knot insertion / removal
- ⌛ Degree elevation / reduction
- ⌛ Construction of common geometries

Basically, chapters 2-4 of [[1]](@ref refs) are implemented so far (adapted to 1-based indexing).

!!! note
    Open knot vectors are assumed everywhere, if not stated otherwise.

!!! note
    The parametric space is assumed to be ``[0, 1]`` or ``{[0,1]}^2``  everywhere.


---
## Installation

Installing NURBS is done by entering the package manager (enter `]` at the julia REPL) and issuing:

```
pkg> add NURBS 
```

---
## [References](@id refs)

The implementation is based on
- [1] L. Piegl, *The NURBS Book*, Berlin Heidelberg, Springer-Verlag, 1997.
- [2] R.N. Simpson, et. al, *A Two-Dimensional Isogeometric Boundary Element Method for Elastostatic Analysis*, Comput. Methods Appl. Mech. Engrg., 2012.
