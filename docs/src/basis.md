
# Bases

The considered 
- B-spline, 
- Curry-Schoenberg, and 
- NURBS 
basis functions are defined by initializing corresponding structures. A naive evaluation method is available for all, an efficient evaluation only for B-splines and Curry-Schoenberg splines.

---
## Define Structures

B-spline, Curry-Schoenberg, and NURBS bases are defined by initializing a [`Bspline`](@ref Bspline), a [`CurrySchoenberg`](@ref CurrySchoenberg), or a [`NURB`](@ref NURB) structure, respectively.

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
Cspl = CurrySchoenberg(p, kVec)
nothing # hide
```


---
## Naive Evaluation

To evaluate the bases at certain points the [`evalNaive`](@ref evalNaive) function is provided. 
It directly implemens the defining equations of the [B-splines](@ref bspl), the [CurrySchoenberg](@ref csspl), and [NURBS](@ref nurbs).
For the derivatives the [`evalNaiveDerivative`](@ref evalNaiveDerivative) function is provided derectly implementing the defining equations of the [derivatives](@ref derB).

!!! note
    The naive evaluation methods are solely implemented to play around with parameters (to get familiar with NURBS and B-splines). 

```@example basis
# --- define points where the bases are evaluated
evalpoints = collect(0:0.001:1.0)

# --- evaluate bases (4-th basis function)
bspline = evalNaive(Bspl, 4, evalpoints) 
nurb    = evalNaive(Nrbs, 4, evalpoints)
cspl    = evalNaive(Cspl, 4, evalpoints)

# --- evaluate derivatives (1st derivative of 4-th basis function)
bsplineD = evalNaiveDerivative(Bspl, 4, 1, evalpoints) 
nurbD    = evalNaiveDerivative(Nrbs, 4, 1, evalpoints)


using Plots
plotly()

Plots.plot(evalpoints, bspline, w=2, 
    label="B-spline", 
    title="4-th basis function", 
    xlabel="𝑢", 
    ylabel="𝑁₄,₂(𝑢)")
Plots.plot!(evalpoints, nurb, w=2, label="NURB")
Plots.plot!(evalpoints, cspl, w=2, label="Curry-Schoenberg")
xlims!(0, 1) # hide
ylims!(0, 1) # hide
savefig("plotBspl.html"); nothing # hide
```

```@raw html
<object data="../plotBspl.html" type="text/html"  style="width:120%;height:50vh;"> </object>
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

For the B-splines and the Curry-Schoenberg splines the efficient evaluation of [[1]](@ref refs) is implemented by handing evaluation points to the [`Bspline`](@ref Bspline) structure and optionally as second argument (only for the B-splines) the maximum derivative that shall be computed. That is, only the basis functions different from zero are evaluated:

!!! note
    For the evaluation of NURBS curves and surfaces (and their derivatives) the B-spline evaluation is sufficient.

```@example basis
bspline  = Bspl(evalpoints)
crySbrg  = Cspl(evalpoints)
bsplineD = Bspl(evalpoints, 2) # 0th, 1st, and 2nd derivative


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
<object data="../plotBspleff.html" type="text/html"  style="width:120%;height:50vh;"> </object>
```

```@example basis
Plots.plot(evalpoints, crySbrg, w=2, 
    leg=false, 
    title="all basis functions", 
    xlabel="𝑢", 
    ylabel="𝑛ᵢ,₂(𝑢)")
xlims!(0, 1) # hide
ylims!(0, 2.0) # hide
savefig("plotCspleff.html"); nothing # hide
```

```@raw html
<object data="../plotCspleff.html" type="text/html"  style="width:120%;height:50vh;"> </object>
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


### Memory Preallocation

In case B-splines shall be evaluated repeatedly memory can be preallocated and reused in subsequent calls:

!!! tip
    This is useful especially if single points are evaluated many times on demand.

```@example basis
pM = NURBS.preAlloc(p, evalpoints)
pD = NURBS.preAllocDer(p, evalpoints, 2)

bspline  = Bspl(evalpoints, pM)
crySbrg  = Cspl(evalpoints, pM)
bsplineD = Bspl(evalpoints, 2, pD) # 0th, 1st, and 2nd derivative
nothing # hide
```
