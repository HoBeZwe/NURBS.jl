
"""
    nurbsNaiveDerivative(basis::NURB, i::Int, k::Int, evalpoints)

Compute the k-th derivative of i-th NURBS basis function evaluated at all 'evalpoints'.
"""
function nurbsNaiveDerivative(basis::NURB, i::Int, k::Int, evalpoints)

    k < 1 && error("The k-th derivative has to be k ≥ 1!")
    k > basis.degree && return zeros(size(evalpoints)) # p+1 th derivative is zero

    return nurbsNaiveDerivative(basis.knotVec, i, basis.degree, basis.weights, evalpoints, k)
end


"""
    nurbsNaiveDerivative(knotVector, i::Int, degree::Int, weights, evalpoints, k::Int; normalize=true)

Compute the k-th derivative of i-th NURBS basis function of degree 'degree' evaluated at all 'evalpoints'.
"""
function nurbsNaiveDerivative(knotVector, i::Int, degree::Int, weights, evalpoints, k::Int)

    # array to store the evaluated points in
    R = similar(evalpoints)

    # loop over all points to be evaluated
    for (ind, u) in enumerate(evalpoints)
        R[ind] = nurbsNaiveDerivative(knotVector, i, degree, weights, u, k)
    end

    return R
end


"""
    nurbsNaiveDerivative(knotVector, i::Int, degree::Int, weights, u::Real, k::Int)

k-th derivative of the i-th NURBS basis function of degree 'degree' evaluated at (single) point 'u'.

Formula (52) of [2].
"""
function nurbsNaiveDerivative(knotVector, i::Int, degree::Int, weights, u::Real, k::Int)

    # handle 0-th derivative case (end of recursion)
    k == 0 && return nurbsNaive(knotVector, i, degree, u, weights)

    # handle higher degree cases
    numBasis = length(knotVector) - degree - 1
    W = 0.0
    for ind in 1:numBasis
        W += bSplineNaive(knotVector, ind, degree, u) * weights[ind]
    end

    R = weights[i] * bSplineNaiveDerivative(knotVector, i, degree, u, k)

    for j in 1:k
        R -=
            binomial(k, j) *
            normalizationDerivative(knotVector, degree, weights, u, j) *
            nurbsNaiveDerivative(knotVector, i, degree, weights, u, k - j)
    end

    return R / W
end


"""
    normalizationDerivative(knotVector, degree::Int, weights, u::Real, k::Int)

Return the k-th derivative of the weight function (21) in [2].
"""
function normalizationDerivative(knotVector, degree::Int, weights, u::Real, k::Int)

    numBasis = length(knotVector) - degree - 1
    W = 0.0
    for ind in 1:numBasis
        W += bSplineNaiveDerivative(knotVector, ind, degree, u, k) * weights[ind]
    end

    return W
end
