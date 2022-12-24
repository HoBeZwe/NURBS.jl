# NURBS

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://hobezwe.github.io/NURBS.jl/dev/)
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/hobezwe/NURBS.jl/blob/main/LICENSE)
[![Build Status](https://github.com/hobezwe/NURBS.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/hobezwe/NURBS.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/HoBeZwe/NURBS.jl/branch/main/graph/badge.svg?token=4F9NUNRC1K)](https://codecov.io/gh/HoBeZwe/NURBS.jl)


## Introduction

This package provides B-spline and NURBS basis functions, their derivatives, as well as curves and surfaces based on both considered basis functions.


## Crucial Assumptions

Open knot vectors are assumed everywhere, if not stated otherwise.

The parametric space is assumed to be [0, 1] everywhere.


## References

The implementation is based on
- L. Piegl, *The NURBS Book*, Berlin Heidelberg, Springer-Verlag, 1997.
- R.N. Simpson, et. al, *A Two-Dimensional Isogeometric Boundary Element Method for Elastostatic Analysis*, Comput. Methods Appl. Mech. Engrg., 2012.


## Documentation

Here you can find the [documentation](https://hobezwe.github.io/NURBS.jl/dev/).