
# Bases

The considered B-spline and NURBS basis functions are defined by initializing corresponding structures. A naive evaluation method is available for both, an efficient evaluation only for B-splines.

---
## Define Structures

NURBS and B-spline bases are defined by initializing a [`Bspline`](@ref Bspline) or a [`NURB`](@ref NURB) structure, respectively.

```@example basis
using NURBS # hide
# --- parameters
b = 6       # number of basis functions
p = 2       # degree of NURBS

w = ones(b) # weights for NURBS
w[4] = 1.8

# --- resulting knot vector
kVec = generateKnotVec(b, p)

# --- structures
Bspl = Bspline(p, kVec)
Nrbs = NURB(p, kVec, w)
nothing # hide
```


---
## Naive Evaluation

To evaluate the bases at certain points the [`bSplineNaive`](@ref bSplineNaive) and the [`nurbsNaive`](@ref nurbsNaive) functions are provided. 
These directly implement the defining equations of the [B-splines](@ref bspl) and [NURBS](@ref nurbs).
For the derivatives the [`bSplineNaiveDerivative`](@ref bSplineNaive) and the [`nurbsNaiveDerivative`](@ref bSplineNaive) functions are provided derectly implementing the defining equations of the [derivatives](@ref derB).

!!! note
    The naive evaluation methods are solely implemented to play around with parameters (to get familiar with NURBS and B-splines). 

```@example basis
# --- define points where the bases are evaluated
evalpoints = collect(0:0.001:1.0)

# --- evaluate bases (4-th basis function)
bspline = bSplineNaive(Bspl, 4, evalpoints) 
nurb    = nurbsNaive(Nrbs, 4, evalpoints)

# --- evaluate derivatives (1st derivative of 4-th basis function)
bsplineD = bSplineNaiveDerivative(Bspl, 4, 1, evalpoints) 
nurbD    = nurbsNaiveDerivative(Nrbs, 4, 1, evalpoints)


using Plots
plotly()

Plots.plot(evalpoints, bspline, w=2, 
    label="B-spline", 
    title="4-th basis function", 
    xlabel="𝑢", 
    ylabel="𝑁₄,₂(𝑢)")
Plots.plot!(evalpoints, nurb, w=2, label="NURB")
xlims!(0, 1) # hide
ylims!(0, 1) # hide
savefig("plotBspl.html"); nothing # hide
```

```@raw html
<object data="plotBspl.html" type="text/html"  style="width:120%;height:50vh;"> </object>
```


```@example basis
Plots.plot(evalpoints, bsplineD, w=2, 
    label="B-spline", 
    title="1st derivative of 4-th basis function", 
    xlabel="𝑢", 
    ylabel="∂ᵤ 𝑁₄,₂(𝑢)")
Plots.plot!(evalpoints, nurbD, w=2, label="NURB")
xlims!(0, 1) # hide
savefig("plotBsplD.html"); nothing # hide
```

```@raw html
<object data="../plotBsplD.html" type="text/html"  style="width:120%;height:50vh;"> </object>
```


---
## Efficient Evaluation

For the B-splines the efficient evaluation of [[1]](@ref refs) is implemented via [`bSpline`](@ref bSpline) and [`bSplineDerivatives`](@ref bSplineDerivatives). That is, only the basis functions different from zero are evaluated:

!!! note
    For the evaluation of NURBS curves and surfaces (and their derivatives) the B-spline evaluation is sufficient.

```@example basis
bspline  = bSpline(Bspl, evalpoints)
bsplineD = bSplineDerivatives(Bspl, 1, evalpoints)

Plots.plot(evalpoints, bspline, w=2, 
    leg=false, 
    title="all basis functions", 
    xlabel="𝑢", 
    ylabel="𝑁ᵢ,₂(𝑢)")
xlims!(0, 1) # hide
ylims!(0, 1) # hide
savefig("plotBspleff.html"); nothing # hide
```

```@raw html
<object data="plotBspleff.html" type="text/html"  style="width:120%;height:50vh;"> </object>
```

```@example basis
Plots.plot(evalpoints, bsplineD[:, 2, :], w=2, 
    leg=false,
    title="1st derivative of all basis functions", 
    xlabel="𝑢", 
    ylabel="∂ᵤ 𝑁ᵢ,₂(𝑢)")
xlims!(0, 1) # hide
savefig("plotBsplDeff.html"); nothing # hide
```

```@raw html
<object data="../plotBsplDeff.html" type="text/html"  style="width:120%;height:50vh;"> </object>
```