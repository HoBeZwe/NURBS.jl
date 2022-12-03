
"""
    bsplineNaive(knotVector, i::Int, degree::Int, evalpoints; normalize=true)

i-th b-spline basis function of degree 'degree' evaluated at all 'evalpoints'.
"""
function bsplineNaive(knotVector, i::Int, degree::Int, evalpoints; normalize=false)

    # normalize knot vector entries to [0, 1]
    if normalize
        knotVector ./= maximum(knotVector)
        @info "The knot vector is being modified."
    end

    # array to store the evaluated points in
    N = similar(evalpoints)

    # loop over all points to be evaluated
    for (ind, u) in enumerate(evalpoints)
        N[ind] = bsplineNaive(knotVector, i, degree, u)
    end

    return N
end


"""
    bsplineNaive(knotVector, i::Int, degree::Int, u::Real)

i-th b-spline basis function of degree 'degree' evaluated at u.

Formula (2.5) of 'The NURBS Book' p. 50.
"""
function bsplineNaive(knotVector, i::Int, degree::Int, u::Real)

    # handle degree 0 case (end of recursion)
    degree == 0 && return bsplineNaive(knotVector, i, u)

    # handle higher degree cases
    ui   = knotVector[i]
    ui1  = knotVector[i + 1]
    uiP  = knotVector[i + degree]
    uiP1 = knotVector[i + degree + 1]

    coeff1 = (u - ui) / (uiP - ui)
    coeff2 = (uiP1 - u) / (uiP1 - ui1)

    # simple fix for division by 0 (not verified that it works in all scenarios)
    isfinite(coeff1) || (coeff1 = 0.0)
    isfinite(coeff2) || (coeff2 = 0.0)

    N = coeff1 * bsplineNaive(knotVector, i, degree - 1, u) + coeff2 * bsplineNaive(knotVector, i + 1, degree - 1, u)

    return N
end


"""
    bsplineNaive(knotVector, i::Int, u::Real)

i-th b-spline basis function of degree 0 evaluated at u.

Formula (2.5) of 'The NURBS Book' p. 50.
"""
function bsplineNaive(knotVector, i::Int, u::Real)

    ui  = knotVector[i]
    ui1 = knotVector[i + 1]

    if u ≥ ui && u < ui1
        return 1.0
    else
        return 0.0
    end
end




"""
    findspan(n::Int, u, knotVector)

Find the spans of a rational B-spline knot vector at the parametric points 'u', where 'b' is the number of basis functions (control points).

Span: the intervall index in which a point lies. E.g., knotVector = [0, 1, 2, 3, 4]. Hence, there are 4 intervalls. u=1.2 lies in the second intervall.

Modification of Algorithm A2.1 from 'The NURBS Book' p. 68.

Assumption that the knotVector is open! (the first and last knot are repeated degree + 1 times)
"""
function findspan(b::Int, u, knotVector)

    # check input
    (minimum(u) < knotVector[1]) && (maximum(u) > knotVector[end]) && error("Some value is outside the knot span")

    # initialize the vector containing the indices 
    spanVec = similar(u, eltype(b))

    for (j, uEval) in enumerate(u)

        # special case: uEval == to last knot vector entry
        if uEval == knotVector[b + 2]
            spanVec[j] = b
            continue
        end

        spanVec[j] = findlast(knotVector .≤ uEval)
    end

    return spanVec
end


"""
    basisfun(knotSpan, uVector, degree::Int, knotVector)

Compute the nonvanishing basis functions of degree 'degree' at the parametric points defined by 'uVector'

Return the basis functions vector of size length(uVector) * (degree + 1).

Adapted from Algorithm A2.2 from 'The NURBS Book' p. 70.
"""
function basisfun(knotSpan, uVector, degree::Int, knotVector)

    B = zeros(length(uVector), degree + 1)
    N = zeros(degree + 1)

    for (jj, u) in enumerate(uVector)

        i = knotSpan[jj] #+ 1 #findspan uses 0-based numbering

        left  = zeros(degree + 1)
        right = zeros(degree + 1)
        N[1]  = 1

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


"""
    curvePoints(nbasisFun::Int, degree::Int, knotVector, controlPoints, uVector)

Compute a 1D curve: given the 'knotVector', the 'controlPoints', and the 'degree', the curve is evaluated at the points given in 'uVector'.

Example for the controlPoints:
    
P1 = SVector(0.0, 0.0, 0.0)
P2 = SVector(0.1, 0.25, 0.0)
P3 = SVector(0.25, 0.3, 0.0)

controlPoints = [P1, P2, P3]
"""
function curvePoints(degree::Int, knotVector, controlPoints, uVector)

    # the number of basis functions is determined by the number of knot vectors and the degree
    nbasisFun = length(knotVector) - degree - 1

    # determine the basis functions evaluated at uVector
    spans = findspan(nbasisFun, uVector, knotVector)
    N = basisfun(spans, uVector, degree, knotVector)

    # determine the curve values
    curve = [SVector(0.0, 0.0, 0.0) for i in eachindex(uVector)] # initialize

    for (j, span) in enumerate(spans)
        for ind in 1:(degree + 1)

            curve[j] += N[j, ind] * controlPoints[span - degree + ind - 1]
        end
    end

    return curve
end


"""
    surfacePoints(uDegree::Int, vDegree::Int, uKnotVector, vKnotVector, controlPoints, uVector, vVector)

Compute surface: given the knotvectors and the degrees in 'u' and 'v' direction, the surface is evaluated at the evaluation points (uVector, vVector).

Control points ordering P_(xi,yj):

P_11 ----- P_12 ----- P_13 ---> y / v direction
  |          |         |
  |          |         |
P_21 ----- P_22 ----- P_23
  |          |         |
  |          |         |
P_31 ----- P_32 ----- P_33
  |
  x / u direction

"""
function surfacePoints(uDegree::Int, vDegree::Int, uKnotVector, vKnotVector, controlPoints, uVector, vVector)

    # u-direction: determine the basis functions evaluated at uVector 
    nbasisFun = length(uKnotVector) - uDegree - 1
    uSpan = findspan(nbasisFun, uVector, uKnotVector)
    Nu = basisfun(uSpan, uVector, uDegree, uKnotVector)

    # v-direction: determine the basis functions evaluated at vVector
    nbasisFun = length(vKnotVector) - vDegree - 1
    vSpan = findspan(nbasisFun, vVector, vKnotVector)
    Nv = basisfun(vSpan, vVector, vDegree, vKnotVector)

    # intialize
    surface = [SVector(0.0, 0.0, 0.0) for i in eachindex(uVector), j in eachindex(vVector)]

    # determine the surface values
    for uPointInd in eachindex(uVector)

        uind = uSpan[uPointInd] - uDegree - 1

        for vPointInd in eachindex(vVector)

            for i in 1:(vDegree + 1)

                temp = SVector(0.0, 0.0, 0.0)

                vind = vSpan[vPointInd] - vDegree + i - 1
                for k in 1:(uDegree + 1)
                    temp = temp + Nu[uPointInd, k] * controlPoints[uind + k, vind]
                end

                surface[uPointInd, vPointInd] += Nv[vPointInd, i] * temp
            end
        end
    end

    return surface
end
