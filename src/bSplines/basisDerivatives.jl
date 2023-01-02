
"""
    bSplineNaiveDerivative(basis::Bspline, i::Int, k::Int, evalpoints)

Compute the k-th derivative of i-th b-spline basis function evaluated at all 'evalpoints'.
"""
function bSplineNaiveDerivative(basis::Bspline, i::Int, k::Int, evalpoints)

    k < 1 && error("The k-th derivative has to be k ≥ 1!")
    k > basis.degree && return zeros(size(evalpoints)) # p+1 th derivative of a polynomial of degree p is zero

    # normalize knot vector entries to [0, 1]
    minimum(basis.knotVec) != 0.0 && error("The knot vector has to start at 0.")
    if maximum(basis.knotVec) != 1.0
        basis.knotVec ./= maximum(basis.knotVec)
        @info "The knot vector is being modified (normalized)."
    end

    return bSplineNaiveDerivative(basis.knotVec, i, basis.degree, evalpoints, k)
end


"""
    bSplineNaiveDer(knotVector, i::Int, degree::Int, evalpoints, k::Int; normalize=true)

Compute the k-th derivative of i-th b-spline basis function of degree 'degree' evaluated at all 'evalpoints'.
"""
function bSplineNaiveDerivative(knotVector, i::Int, degree::Int, evalpoints, k::Int)

    # array to store the evaluated points in
    N = similar(evalpoints)

    # loop over all points to be evaluated
    for (ind, u) in enumerate(evalpoints)
        N[ind] = bSplineNaiveDerivative(knotVector, i, degree, u, k)
    end

    return N
end


"""
    bSplineNaiveDerivative(knotVector, i::Int, degree::Int, u::Real, k::Int)

k-th derivative of the i-th b-spline basis function of degree 'degree' evaluated at (single) point 'u'.

Formula (2.9) of 'The NURBS Book' p. 61. (Recursive implementation to avoid the faculties in the (2.10) formula.)
"""
function bSplineNaiveDerivative(knotVector, i::Int, degree::Int, u::Real, k::Int)

    # handle 0-th derivative case (end of recursion)
    k == 0 && return bSplineNaive(knotVector, i, degree, u)

    # handle higher degree cases
    ui   = knotVector[i]
    ui1  = knotVector[i + 1]
    uiP  = knotVector[i + degree]
    uiP1 = knotVector[i + degree + 1]

    coeff1 = degree / (uiP - ui)
    coeff2 = degree / (uiP1 - ui1)

    # handle division by 0 ('The NURBS Book' remarks on p. 61)
    isfinite(coeff1) || (coeff1 = 0.0)
    isfinite(coeff2) || (coeff2 = 0.0)

    N =
        coeff1 * bSplineNaiveDerivative(knotVector, i, degree - 1, u, k - 1) -
        coeff2 * bSplineNaiveDerivative(knotVector, i + 1, degree - 1, u, k - 1)

    return N
end


"""
    bSplineDerivatives(basis::Bspline, k::Int, evalpoints)

Evaluate k-the derivative of B-spline basis at all evalpoints (all basis functions different from 0 at the evalpoints are evaluated).
"""
function bSplineDerivatives(basis::Bspline, k::Int, evalpoints)

    k < 0 && error("The k-th derivative has to be k ≥ 0!")

    numBasis = numBasisFunctions(basis)
    knotSpan = findSpan(numBasis, evalpoints, basis.knotVec)

    return derBasisFun(knotSpan, basis.degree, evalpoints, basis.knotVec, k)
end


"""
    derBasisFun(knotSpan, degree::Int, evalpoints, knotVector, numberDerivatives::Int)

Compute the nonvanishing B-spline basis functions and its derivatives of degree 'degree' at the parametric points defined by 'uVector'.

Organization of output:  dersv[n, k, :] contains (k-1)-th derivative at n-th point.
   
Adapted from Algorithm A2.3 from 'The NURBS BOOK' p. 72.
"""
function derBasisFun(knotSpan, degree::Int, evalpoints, knotVector, numberDerivatives::Int)

    dersv = zeros(length(evalpoints), numberDerivatives + 1, degree + 1)

    for jj in eachindex(evalpoints)

        i = knotSpan[jj] #+ 1 # convert to base-1 numbering of knot spans
        u = evalpoints[jj]

        ders  = zeros(numberDerivatives + 1, degree + 1)
        ndu   = zeros(degree + 1, degree + 1)
        left  = zeros(degree + 1)
        right = zeros(degree + 1)
        a     = zeros(2, degree + 1)

        ndu[1, 1] = 1

        for j in 1:degree
            left[j + 1] = u - knotVector[i + 1 - j]
            right[j + 1] = knotVector[i + j] - u
            saved = 0
            for r in 0:(j - 1)
                ndu[j + 1, r + 1] = right[r + 2] + left[j - r + 1]
                temp = ndu[r + 1, j] / ndu[j + 1, r + 1]
                ndu[r + 1, j + 1] = saved + right[r + 2] * temp
                saved = left[j - r + 1] * temp
            end
            ndu[j + 1, j + 1] = saved
        end

        for j in 0:degree
            ders[1, j + 1] = ndu[j + 1, degree + 1]
        end

        for r in 0:degree
            s1 = 0
            s2 = 1
            a[1, 1] = 1
            for k in 1:numberDerivatives # compute kth derivative
                d = 0.0
                rk = r - k
                pk = degree - k

                if r ≥ k
                    a[s2 + 1, 1] = a[s1 + 1, 1] / ndu[pk + 2, rk + 1]
                    d = a[s2 + 1, 1] * ndu[rk + 1, pk + 1]
                end

                if rk ≥ -1
                    j1 = 1
                else
                    j1 = -rk
                end
                if (r - 1) <= pk
                    j2 = k - 1
                else
                    j2 = degree - r
                end

                for j in j1:j2
                    a[s2 + 1, j + 1] = (a[s1 + 1, j + 1] - a[s1 + 1, j]) / ndu[pk + 2, rk + j + 1]
                    d = d + a[s2 + 1, j + 1] * ndu[rk + j + 1, pk + 1]
                end

                if r ≤ pk
                    a[s2 + 1, k + 1] = -a[s1 + 1, k] / ndu[pk + 2, r + 1]
                    d = d + a[s2 + 1, k + 1] * ndu[r + 1, pk + 1]
                end

                ders[k + 1, r + 1] = d
                j = s1
                s1 = s2
                s2 = j
            end
        end

        r = degree
        for k in 1:numberDerivatives
            for j in 0:degree
                ders[k + 1, j + 1] = ders[k + 1, j + 1] * r
            end
            r = r * (degree - k)
        end

        dersv[jj, :, :] = ders
    end

    return dersv
end
