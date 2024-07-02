
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/src/assets/logo_Scat_READMEwhite.svg" height="190">
  <source media="(prefers-color-scheme: light)" srcset="docs/src/assets/logo_Scat_README.svg" height="190">
  <img alt="" src="" height="190">
</picture>

[![Docs-stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://hobezwe.github.io/NURBS.jl/stable/)
[![Docs-dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://hobezwe.github.io/NURBS.jl/dev/)
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/hobezwe/NURBS.jl/blob/main/LICENSE)
[![Build Status](https://github.com/hobezwe/NURBS.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/hobezwe/NURBS.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/HoBeZwe/NURBS.jl/branch/main/graph/badge.svg?token=4F9NUNRC1K)](https://codecov.io/gh/HoBeZwe/NURBS.jl)
[![DOI](https://zenodo.org/badge/579998043.svg)](https://zenodo.org/badge/latestdoi/579998043)



## Introduction

This package provides functionality to define and evaluate B-spline, Curry-Schoenberg, and NURBS (non-uniform rational B-spline) basis functions, their derivatives, as well as curves and surfaces based on B-spline and NURBS basis functions.

The following aspects are implemented (✓) and planned (⌛):

##### B-spline, Curry-Schoenberg & NURBS evaluation
- ✓ Basis & derivatives
- ✓ Curves & derivatives
- ✓ Surfaces & derivatives

##### Fundamental operations
- ✓ File I/O (.step)
- ✓ Knot manipulation
    - knot insertion / refinement
    - knot removal
    - splitting of curves and surfaces
- ✓ Transformation of curves and surfaces
    - scaling
    - translating
    - rotating
    - mirroring
- ⌛ Degree elevation / reduction
- ⌛ Construction of common geometries

##### Connectivity
- ✓ Determine patch connectivity
    - identify interfaces between patches
    - introduce per patch local numbering for vertices and edges
- ✓ Virtual Bezier mesh connectivty (for FEM)
    - introduce on each patch a virtual Bezier mesh
    - determine adjacency information of mesh cells


## Citation

Please cite this package following the information on [Zenodo](https://zenodo.org/badge/latestdoi/579998043).



## Documentation

- Documentation for the [latest stable version](https://hobezwe.github.io/NURBS.jl/stable/).
- Documentation for the [development version](https://hobezwe.github.io/NURBS.jl/dev/).