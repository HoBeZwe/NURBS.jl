

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/src/assets/logo_Scat_READMEwhite.svg" height="190">
  <source media="(prefers-color-scheme: light)" srcset="docs/src/assets/logo_Scat_README.svg" height="190">
  <img alt="" src="" height="190">
</picture>


[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://hobezwe.github.io/Nurbs.jl/dev/)
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/hobezwe/Nurbs.jl/blob/main/LICENSE)
[![Build Status](https://github.com/hobezwe/Nurbs.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/hobezwe/Nurbs.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/HoBeZwe/Nurbs.jl/branch/main/graph/badge.svg?token=4F9NUNRC1K)](https://codecov.io/gh/HoBeZwe/Nurbs.jl)
[![DOI](https://zenodo.org/badge/375493054.svg)](https://zenodo.org/badge/latestdoi/375493054)


## Introduction

This package provides functionality to define and evaluate B-spline and NURBS (non-uniform rational B-spline) basis functions, their derivatives, as well as curves and surfaces based on both considered basis functions.

The following aspects are implemented (✓) and planned (⌛):

##### B-spline & NURBS evaluation
- ✓ Basis & derivatives 
- ✓ Curves & derivatives 
- ✓ Surfaces & derivatives 

##### Fundamental operations
- ✓ File I/O (basic)
- ⌛ Knot insertion / removal
- ⌛ Degree elevation / reduction
- ⌛ Construction of common geometries

## Citation

Please cite this package following the information on [Zenodo](https://zenodo.org/badge/latestdoi/375493054).


## Documentation

Here you can find the [documentation](https://hobezwe.github.io/Nurbs.jl/dev/).