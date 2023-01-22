
"""
    nurbsNaive(basis::NURB, i::Int, evalpoints)

i-th NURB basis function evaluated at all 'evalpoints'.
"""
function nurbsNaive(basis::NURB, i::Int, evalpoints)

    return nurbsNaive(basis.knotVec, i, basis.degree, evalpoints, basis.weights)
end


"""
    nurbsNaive(knotVector, i::Int, degree::Int, evalpoints, weights)

i-th NURB basis function of degree 'degree', the 'knotVector', and with 'weights' evaluated at all 'evalpoints'.
"""
function nurbsNaive(knotVector, i::Int, degree::Int, evalpoints, weights)

    # array to store the evaluated points in
    N = similar(evalpoints)

    # loop over all points to be evaluated
    for (ind, u) in enumerate(evalpoints)
        N[ind] = nurbsNaive(knotVector, i, degree, u, weights)
    end

    return N
end


"""
    nurbsNaive(knotVector, i::Int, degree::Int, u::Real, weights)

i-th NURB basis function of degree 'degree', the 'knotVector', and with 'weights' evaluated at single point 'u'.

Formula (4.2) of 'The NURBS Book' p. 118. 
"""
function nurbsNaive(knotVector, i::Int, degree::Int, u::Real, weights)

    N = bSplineNaive(knotVector, i, degree, u) * weights[i]

    numBasis = length(knotVector) - degree - 1
    normalize = 0.0
    for ind in 1:numBasis
        normalize += bSplineNaive(knotVector, ind, degree, u) * weights[ind]
    end

    return N / normalize
end
