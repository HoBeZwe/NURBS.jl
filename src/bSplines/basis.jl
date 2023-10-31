
"""
    evalNaive(basis::Basis, i::Int, evalpoints)

i-th B-spline or Curry-Schoenberg basis function evaluated at all 'evalpoints'.
"""
function evalNaive(basis::Basis, i::Int, evalpoints)

    norma = normalization(basis, i)

    return norma * bSplineNaive(basis.knotVec, i, basis.degree, evalpoints)
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
    normalization(basis::Basis, i::Int)

For Bsplines there is no normalization.
"""
function normalization(basis::Basis, i::Int)

    return 1.0
end


"""
    normalization(basis::CurrySchoenberg, i::Int)

Normalization for the standard B-splines to obtain Curry-Schoenberg splines.
"""
function normalization(basis::CurrySchoenberg, i::Int)

    return (basis.degree + 1) / (basis.knotVec[i + basis.degree + 1] - basis.knotVec[i])
end


"""
    (basis::Bspline)(evalpoints)

Evaluate B-spline basis at all evalpoints.
"""
function (basis::Union{Bspline,CurrySchoenberg})(evalpoints)

    numBasis = numBasisFunctions(basis)
    knotSpan = findSpan(numBasis, evalpoints, basis.knotVec, basis.degree)

    return basisFun(knotSpan, evalpoints, basis)
end


struct pAlloc{T<:Real,F<:Int}
    B::Matrix{T}
    N::Vector{T}
    left::Vector{T}
    right::Vector{T}
    spanVec::Vector{F}
end


"""
    (basis::Bspline)(evalpoints)

Evaluate B-spline basis at all evalpoints.
"""
function (basis::Union{Bspline,CurrySchoenberg})(evalpoints, prealloc::pAlloc)

    numBasis = numBasisFunctions(basis)
    knotSpan = findSpan!(prealloc.spanVec, numBasis, evalpoints, basis.knotVec, basis.degree)

    return basisFun!(prealloc, knotSpan, evalpoints, basis)
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

    spanVec = similar(uVector, eltype(degree))

    return pAlloc(B, N, left, right, spanVec)
end


"""
    findSpan(b::T, uVec, kVec, p::T)

Allocate memory and call findSpan!    
"""
function findSpan(b::T, uVec, kVec, p::T) where {T<:Int}

    spanVec = similar(uVec, eltype(b))

    return findSpan!(spanVec, b, uVec, kVec, p)
end


"""
    findSpan!(spanVec, b::T, uVec, kVec, p::T)

Find the spans of a B-spline knot vector at the parametric points 'u', where 'b' is the number of basis functions (control points).

Span: the intervall index in which a point lies. E.g., knotVector = [0, 1, 2, 3, 4]. Hence, there are 4 intervalls. u=1.2 lies in the second intervall.

Modification of Algorithm A2.1 from 'The NURBS Book' p. 68.

Assumption that the knotVector is open! (the first and last knot are repeated degree + 1 times)
"""
function findSpan!(spanVec, b::T, uVec, kVec, p::T) where {T<:Int}

    for (j, u) in enumerate(uVec)

        # --- Special case 
        if u == kVec[end]
            spanVec[j] = b
            continue
        end

        # --- Do binary search 
        low  = p + 1
        high = b + 1

        mid = Base.midpoint(low, high)

        while u < kVec[mid] || u >= kVec[mid + 1] # is u in interval [ kVec[mid]...kVec[mid+1] ) ?

            if u < kVec[mid]
                high = mid
            else
                low = mid
            end

            mid = Base.midpoint(low, high)
        end

        spanVec[j] = mid
    end

    return spanVec
end


"""
    basisFun(knotSpan, uVector, basis::Basis)

Allocate memory and call basisFun!    
"""
function basisFun(knotSpan, uVector, basis::Basis)

    prealloc = preAlloc(basis.degree, uVector)

    return basisFun!(prealloc, knotSpan, uVector, basis)
end


"""
    basisFun!(prealloc::pAlloc, knotSpan, uVector, basis::Basis)

Compute the nonvanishing B-spline basis functions of degree 'degree' at the parametric points defined by 'uVector'

Return the B-spline basis functions vector of size length(uVector) * (degree + 1).

Adapted from Algorithm A2.2 from 'The NURBS Book' p. 70.
"""
function basisFun!(prealloc::pAlloc, knotSpan, uVector, basis::Basis)

    degree = basis.degree
    knotVector = basis.knotVec

    B = prealloc.B
    N = prealloc.N

    left  = prealloc.left
    right = prealloc.right

    for (jj, u) in enumerate(uVector)

        i = knotSpan[jj] #+ 1 #findspan uses 0-based numbering
        N[1] = 1.0

        for j in 1:degree
            left[j + 1]  = u - knotVector[i + 1 - j]
            right[j + 1] = knotVector[i + j] - u
            saved        = 0.0

            for r in 1:j
                temp = N[r] / (right[r + 1] + left[j - r + 2])
                N[r] = saved + right[r + 1] * temp
                saved = left[j - r + 2] * temp
            end
            N[j + 1] = saved
        end

        normalize!(N, i, degree, knotVector, basis)

        B[jj, :] = N
    end

    return B
end


"""
    normalize!(N, i::T, degree::T, knotVector, basis::Bspline)

For Bsplines there is no normalization.
"""
function normalize!(N, i::T, degree::T, knotVector, basis::Basis) where {T<:Int}

    return nothing
end


"""
    normalize!(N, i::T, degree::T, knotVector, basis::CurrySchoenberg)

Normalize the standard B-splines to obtain Curry-Schoenberg splines.
"""
function normalize!(N, i::T, degree::T, knotVector, basis::CurrySchoenberg) where {T<:Int}

    dp1 = degree + 1

    for j in 1:dp1
        ind = i - dp1 + j # index follows from the knotspan
        N[j] *= dp1 / (knotVector[ind + dp1] - knotVector[ind])
    end
end
