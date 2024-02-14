
# NURBS.jl

This package provides functionality to define and evaluate B-spline, Curry-Schoenberg, and NURBS (non-uniform rational B-spline) basis functions, their derivatives, as well as curves and surfaces based on B-spline and NURBS basis functions.


---
## Overview

The following aspects are implemented (✓) and planned (⌛):

##### B-spline, Curry-Schoenberg & NURBS evaluation
- ✓ Basis & derivatives
- ✓ Curves & derivatives
- ✓ Surfaces & derivatives

##### Fundamental operations
- ✓ File I/O (.step)
- ⌛ Knot manipulation
    - ✓ knot insertion / refinement
    - ⌛ knot removal
    - ✓ splitting of curves and surfaces
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
- [3] C. de Boor, *A Practical Guide to Splines*, revised ed., Appl. Math. Sci., vol. 27, Springer-Verlag, New York, 2001.
- [4] L. Beirão da Veiga, A. Buffa, G. Sangalli, R. Vázquez, *Analysis-suitable T-splines of arbitrary degree: Definition, linear independence and approximation properties*, Math. Models Methods Appl. Sci. 23, 2013.
