
"""
    evalNaive(basis::Basis, i::Int, evalpoints)

i-th B-spline or Curry-Schoenberg basis function evaluated at all 'evalpoints'.
"""
function evalNaive(basis::Basis, i::Int, evalpoints)

    norma = normalization(basis, i)

    return norma * bSplineNaive(basis.knotVec, i, degree(basis), evalpoints)
end


"""
    bSplineNaive(knotVector, i::Int, degree::Int, evalpoints)

i-th b-spline basis function of degree 'degree' evaluated at all 'evalpoints'.

The knotvector is assumed to be normalized to [1, 0].
"""
function bSplineNaive(knotVector, i::Int, degree::Int, evalpoints::AbstractArray{T}) where {T}

    # array to store the evaluated points in
    B = similar(evalpoints)

    N = numBasisFunctions(knotVector, degree)

    # loop over all points to be evaluated
    for (ind, u) in enumerate(evalpoints)
        B[ind] = bSplineNaive(knotVector, i, degree, u, N)
    end

    return B
end


"""
    bSplineNaive(knotVector, i::Int, degree::Int, u::Real, N)

i-th b-spline basis function of degree 'degree' evaluated at u.

Formula (2.5) of 'The NURBS Book' p. 50.
"""
function bSplineNaive(knotVector, i::Int, degree::Int, u::T, N=numBasisFunctions(knotVector, degree)) where {T}

    i == N && u == one(T) && return one(T)

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
    isfinite(coeff1) || (coeff1 = zero(T))
    isfinite(coeff2) || (coeff2 = zero(T))

    B = coeff1 * bSplineNaive(knotVector, i, degree - 1, u, N) + coeff2 * bSplineNaive(knotVector, i + 1, degree - 1, u, N)

    return B
end


"""
    bSplineNaive(knotVector, i::Int, u::Real)

i-th b-spline basis function of degree 0 evaluated at u.

Formula (2.5) of 'The NURBS Book' p. 50.
"""
function bSplineNaive(knotVector, i::Int, u::T) where {T<:Real}

    ui  = knotVector[i]
    ui1 = knotVector[i + 1]

    if u ≥ ui && u < ui1
        return one(T)
    else
        return zero(T)
    end
end


"""
    normalization(basis::Basis, i::Int)

For Bsplines there is no normalization.
"""
function normalization(basis::Basis{F}, i::Int) where {F}

    return one(F)
end


"""
    normalization(basis::CurrySchoenberg, i::Int)

Normalization for the standard B-splines to obtain Curry-Schoenberg splines.
"""
function normalization(basis::CurrySchoenberg, i::Int)

    return (degree(basis) + 1) / (basis.knotVec[i + degree(basis) + 1] - basis.knotVec[i])
end


"""
    (basis::Bspline)(evalpoints)

Evaluate B-spline basis at all evalpoints.
"""
function (basis::Union{Bspline,CurrySchoenberg})(evalpoints)

    prealloc = preAlloc(degree(basis), evalpoints)
    basis(evalpoints, prealloc)

    return prealloc.B
end


struct pAlloc{T,F}
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
    findSpan!(prealloc.spanVec, numBasis, evalpoints, basis.knotVec, degree(basis))
    basisFun!(prealloc, prealloc.spanVec, evalpoints, basis)

    return prealloc.B
end


"""
    preAlloc(degree::Int, uVector::Vector{T})

Allocate memory for basisFun.
"""
function preAlloc(degree, uVector)

    F = eltype(uVector)
    p = degree + 1

    B = zeros(F, length(uVector), p)
    N = zeros(F, p)

    left  = zeros(F, p)
    right = zeros(F, p)

    spanVec = allocSpan(uVector, degree)

    return pAlloc(B, N, left, right, spanVec)
end

allocSpan(uVector::Vector{F}, degree::T) where {T,F} = similar(uVector, T)
allocSpan(u::F, degree::T) where {T,F} = return zeros(T, 1)


"""
    findSpan(b::T, uVec, kVec, p::T)

Allocate memory and call findSpan!    
"""
function findSpan(b::T, uVec::AbstractArray, kVec, p::T) where {T}

    spanVec = similar(uVec, T)
    findSpan!(spanVec, b, uVec, kVec, p)

    return spanVec
end


"""
    findSpan!(spanVec, b::T, uVec, kVec, p::T)

Find the spans of a B-spline knot vector at the parametric points 'u', where 'b' is the number of basis functions (control points).

Span: the intervall index in which a point lies. E.g., knotVector = [0, 1, 2, 3, 4]. Hence, there are 4 intervalls. u=1.2 lies in the second intervall.

Modification of Algorithm A2.1 from 'The NURBS Book' p. 68.

Assumption that the knotVector is open! (the first and last knot are repeated degree + 1 times)
"""
function findSpan!(spanVec::Vector{T}, b::T, uVec, kVec, p::T) where {T}

    for j in eachindex(uVec)
        spanVec[j] = findSpan(b, uVec[j], kVec, p)
    end

    return nothing
end


function findSpan(b::T, u::F, kVec, p::T) where {F,T}

    # --- Special case 
    u == kVec[end] && return b

    # --- Do binary search 
    low  = p + 1
    high = b + 1

    mid = Base.midpoint(low, high)

    while u < kVec[mid] || u >= kVec[mid + 1] # is u outside the interval [ kVec[mid]...kVec[mid+1] ) ?

        if u < kVec[mid]
            high = mid
        else
            low = mid
        end

        mid = Base.midpoint(low, high)
    end # if not, we have found the interval

    return mid
end


"""
    basisFun(knotSpan, uVector, basis::Basis)

Allocate memory and call basisFun!    
"""
function basisFun(knotSpan, uVector, basis::Basis)

    prealloc = preAlloc(degree(basis), uVector)
    basisFun!(prealloc, knotSpan, uVector, basis)

    return prealloc.B
end


"""
    basisFun!(prealloc::pAlloc, knotSpan, uVector, basis::Basis)

Compute the nonvanishing B-spline basis functions of degree 'degree' at the parametric points defined by 'uVector'

Return the B-spline basis functions vector of size length(uVector) * (degree + 1).

Adapted from Algorithm A2.2 from 'The NURBS Book' p. 70.
"""
function basisFun!(prealloc::pAlloc, knotSpan, uVector, basis::Basis)

    degr = degree(basis)
    knotVector = basis.knotVec

    B = prealloc.B
    N = prealloc.N

    left  = prealloc.left
    right = prealloc.right

    for (jj, u) in enumerate(uVector)

        i = knotSpan[jj] #+ 1 #findspan uses 0-based numbering
        N[1] = 1.0

        for j in 1:degr
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

        normalize!(N, i, degr, knotVector, basis)

        B[jj, :] = N
    end

    return nothing
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

    return nothing
end
