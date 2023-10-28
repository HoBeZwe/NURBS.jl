
"""
    evalNaive(basis::Bspline, i::Int, evalpoints)

i-th B-spline basis function evaluated at all 'evalpoints'.
"""
function evalNaive(basis::Bspline, i::Int, evalpoints)

    return bSplineNaive(basis.knotVec, i, basis.degree, evalpoints)
end


"""
    bSplineNaive(knotVector, i::Int, degree::Int, evalpoints; normalize=true)

i-th b-spline basis function of degree 'degree' evaluated at all 'evalpoints'.

The knotvector is assumed to be normalized to [1, 0].
"""
function bSplineNaive(knotVector, i::Int, degree::Int, evalpoints)

    # array to store the evaluated points in
    N = similar(evalpoints)

    # loop over all points to be evaluated
    for (ind, u) in enumerate(evalpoints)
        N[ind] = bSplineNaive(knotVector, i, degree, u)
    end

    return N
end


"""
    bSplineNaive(knotVector, i::Int, degree::Int, u::Real)

i-th b-spline basis function of degree 'degree' evaluated at u.

Formula (2.5) of 'The NURBS Book' p. 50.
"""
function bSplineNaive(knotVector, i::Int, degree::Int, u::Real)

    # handle degree 0 case (end of recursion)
    degree == 0 && return bSplineNaive(knotVector, i, u)

    # handle higher degree cases
    ui   = knotVector[i]
    ui1  = knotVector[i + 1]
    uiP  = knotVector[i + degree]
    uiP1 = knotVector[i + degree + 1]

    coeff1 = (u - ui) / (uiP - ui)
    coeff2 = (uiP1 - u) / (uiP1 - ui1)

    # handle division by 0 ('The NURBS Book' p. 51)
    isfinite(coeff1) || (coeff1 = 0.0)
    isfinite(coeff2) || (coeff2 = 0.0)

    N = coeff1 * bSplineNaive(knotVector, i, degree - 1, u) + coeff2 * bSplineNaive(knotVector, i + 1, degree - 1, u)

    return N
end


"""
    bSplineNaive(knotVector, i::Int, u::Real)

i-th b-spline basis function of degree 0 evaluated at u.

Formula (2.5) of 'The NURBS Book' p. 50.
"""
function bSplineNaive(knotVector, i::Int, u::Real)

    ui  = knotVector[i]
    ui1 = knotVector[i + 1]

    if u ≥ ui && u < ui1
        return 1.0
    else
        return 0.0
    end
end


"""
    (basis::Bspline)(evalpoints)

Evaluate B-spline basis at all evalpoints.
"""
function (basis::Bspline)(evalpoints)

    numBasis = numBasisFunctions(basis)
    knotSpan = findSpan(numBasis, evalpoints, basis.knotVec)

    return basisFun(knotSpan, evalpoints, basis.degree, basis.knotVec)
end


struct pAlloc{T<:Real}
    B::Matrix{T}
    N::Vector{T}
    left::Vector{T}
    right::Vector{T}
end


"""
    (basis::Bspline)(evalpoints)

Evaluate B-spline basis at all evalpoints.
"""
function (basis::Bspline)(evalpoints, prealloc::pAlloc)

    numBasis = numBasisFunctions(basis)
    knotSpan = findSpan(numBasis, evalpoints, basis.knotVec)

    return basisFun!(prealloc, knotSpan, evalpoints, basis.degree, basis.knotVec)
end


"""
    preAlloc(degree::Int, uVector::Vector{T})

Allocate memory for basisFun.
"""
function preAlloc(degree::Int, uVector::Vector{T}) where {T}

    B = zeros(T, length(uVector), degree + 1)
    N = zeros(T, degree + 1)

    left  = zeros(T, degree + 1)
    right = zeros(T, degree + 1)

    return pAlloc(B, N, left, right)
end


"""
    findSpan(n::Int, u, knotVector)

Find the spans of a B-spline knot vector at the parametric points 'u', where 'b' is the number of basis functions (control points).

Span: the intervall index in which a point lies. E.g., knotVector = [0, 1, 2, 3, 4]. Hence, there are 4 intervalls. u=1.2 lies in the second intervall.

Modification of Algorithm A2.1 from 'The NURBS Book' p. 68.

Assumption that the knotVector is open! (the first and last knot are repeated degree + 1 times)
"""
function findSpan(b::Int, u, knotVector)

    # check input
    #(minimum(u) < knotVector[1]) && (maximum(u) > knotVector[end]) && error("Some value is outside the knot span")

    # initialize the vector containing the indices 
    spanVec = similar(u, eltype(b))

    for (j, uEval) in enumerate(u)

        # special case: uEval == to last knot vector entry
        if uEval == knotVector[end]
            spanVec[j] = b
            continue
        end

        spanVec[j] = findlast(knotVector .≤ uEval)
    end

    return spanVec
end


"""
    basisFun(knotSpan, uVector, degree::Int, knotVector)

Allocate memory and call basisFun!    
"""
function basisFun(knotSpan, uVector, degree::Int, knotVector)

    prealloc = preAlloc(degree, uVector)

    return basisFun!(prealloc, knotSpan, uVector, degree, knotVector)
end


"""
    basisFun!(prealloc::pAlloc, knotSpan, uVector, degree::Int, knotVector)

Compute the nonvanishing B-spline basis functions of degree 'degree' at the parametric points defined by 'uVector'

Return the B-spline basis functions vector of size length(uVector) * (degree + 1).

Adapted from Algorithm A2.2 from 'The NURBS Book' p. 70.
"""
function basisFun!(prealloc::pAlloc, knotSpan, uVector, degree::Int, knotVector)

    B = prealloc.B
    N = prealloc.N

    left  = prealloc.left
    right = prealloc.right

    for (jj, u) in enumerate(uVector)

        i = knotSpan[jj] #+ 1 #findspan uses 0-based numbering
        N[1] = 1

        for j in 1:degree
            left[j + 1]  = u - knotVector[i + 1 - j]
            right[j + 1] = knotVector[i + j] - u
            saved        = 0

            for r in 0:(j - 1)
                temp = N[r + 1] / (right[r + 2] + left[j - r + 1])
                N[r + 1] = saved + right[r + 2] * temp
                saved = left[j - r + 1] * temp
            end
            N[j + 1] = saved
        end

        B[jj, :] = N
    end

    return B
end
