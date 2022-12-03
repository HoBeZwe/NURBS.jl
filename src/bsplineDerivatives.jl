
"""
    bsplineNaiveDer(knotVector, i::Int, degree::Int, evalpoints, k::Int; normalize=true)

Compute the k-th derivative of i-th b-spline basis function of degree 'degree' evaluated at all 'evalpoints'.
"""
function bsplineNaiveDerivative(knotVector, i::Int, degree::Int, evalpoints, k::Int; normalize=false)

    k < 1 && error("The k-th derivative has to be k ≥ 1!")  
    k > degree && return zeros(size(evalpoints)) # p+1 th derivative of a polynomial of degree p is zero

    # normalize knot vector entries to [0, 1]
    if normalize
        knotVector ./= maximum(knotVector)
        @info "The knot vector is being modified."
    end

    # array to store the evaluated points in
    N = similar(evalpoints)

    # loop over all points to be evaluated
    for (ind, u) in enumerate(evalpoints)
        N[ind] = bsplineNaiveDerivative(knotVector, i, degree, u, k)
    end

    return N
end


"""
    bsplineNaiveDerivative(knotVector, i::Int, degree::Int, u::Real, k::Int)

k-th derivative of the i-th b-spline basis function of degree 'degree' evaluated at (single) point 'u'.

Formula (2.9) of 'The NURBS Book' p. 61. (Recursive implementaion to avoid the faculties in the (2.10) formula.)
"""
function bsplineNaiveDerivative(knotVector, i::Int, degree::Int, u::Real, k::Int)

    # handle 0-th derivative case (end of recursion)
    k == 0 && return bsplineNaive(knotVector, i, degree, u)

    # handle higher degree cases
    ui   = knotVector[i]
    ui1  = knotVector[i + 1]
    uiP  = knotVector[i + degree]
    uiP1 = knotVector[i + degree + 1]

    coeff1 = degree / (uiP - ui)
    coeff2 = degree / (uiP1 - ui1)

    # simple fix for division by 0 
    isfinite(coeff1) || (coeff1 = 0.0)
    isfinite(coeff2) || (coeff2 = 0.0)

    N = coeff1 * bsplineNaiveDerivative(knotVector, i, degree - 1, u, k - 1) - coeff2 * bsplineNaiveDerivative(knotVector, i + 1, degree - 1, u, k - 1)

    return N
end