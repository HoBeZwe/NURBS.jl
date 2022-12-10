
"""
    numBasisFunctions(basis::Basis)

The number of basis functions is fixed by the knot vector and the degree.

Assumption: the first and last knot vector entry has mulitplicity degree + 1.
"""
numBasisFunctions(basis::Basis) = length(basis.knotVec) - basis.degree - 1


"""
    generateKnotVec(b::Int, degree::Int)

Convenience function to generate a knot vector for 'b' basis functions and a certain 'degree': 

The first and last entry are repeated 'degree'+1 times. Normalized to [0, 1].
"""
function generateKnotVec(b::Int, degree::Int)

    mpc = degree + 1     # required mulitplicity of first and last knot
    m = b + degree + 1   # length of knot vector

    aa = zeros(Int, mpc, 1)
    bb = ones(Int, mpc, 1) * (m - 2 * mpc + 1)

    kVec = Vector{Float64}(undef, m)
    kVec[1:mpc] .= aa
    kVec[(end - mpc + 1):end] .= bb
    kVec[(mpc + 1):(end - mpc)] .= collect(1:(bb[1] - 1))

    kVec ./= maximum(kVec)

    return kVec
end


"""
    bSplineNaive(basis::Bspline, i::Int, evalpoints)

i-th B-spline basis function evaluated at all 'evalpoints'.
"""
function bSplineNaive(basis::Bspline, i::Int, evalpoints)

    # normalize knot vector entries to [0, 1]
    minimum(basis.knotVec) != 0.0 && error("The knot vector has to start at 0.")
    if maximum(basis.knotVec) != 1.0
        knotVector ./= maximum(knotVector)
        @info "The knot vector is being modified (normalized)."
    end

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
    bSpline(basis::Bspline, evalpoints)

Evaluate B-spline basis at all evalpoints.
"""
function bSpline(basis::Bspline, evalpoints)

    numBasis = numBasisFunctions(basis)
    knotSpan = findSpan(numBasis, evalpoints, basis.knotVec)

    return basisFun(knotSpan, evalpoints, basis.degree, basis.knotVec)
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
    basisFun(knotSpan, uVector, degree::Int, knotVector)

Compute the nonvanishing B-spline basis functions of degree 'degree' at the parametric points defined by 'uVector'

Return the B-spline basis functions vector of size length(uVector) * (degree + 1).

Adapted from Algorithm A2.2 from 'The NURBS Book' p. 70.
"""
function basisFun(knotSpan, uVector, degree::Int, knotVector)

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
    curvePoints(curve::BsplineCurve, uVector)

Convenience function to plot a NURBS curve.
"""
curvePoints(curve::BsplineCurve, uVector) = curvePoints(curve.basis.degree, curve.basis.knotVec, curve.controlPoints, uVector)


"""
    curvePoints(nbasisFun::Int, degree::Int, knotVector, controlPoints, uVector)

Compute a 1D B-spline curve: given the 'knotVector', the 'controlPoints', and the 'degree', the curve is evaluated at the points given in 'uVector'.

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
    spans = findSpan(nbasisFun, uVector, knotVector)
    N = basisFun(spans, uVector, degree, knotVector)

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
    surfacePoints(Patch::BsplineSurface, uEvalpoints, vEvalpoints)

Convenience function to plot a B-spline surface.
"""
surfacePoints(Patch::BsplineSurface, uEvalpoints, vEvalpoints) = surfacePoints(
    Patch.uBasis.degree,
    Patch.vBasis.degree,
    Patch.uBasis.knotVec,
    Patch.vBasis.knotVec,
    Patch.controlPoints,
    uEvalpoints,
    vEvalpoints,
)


"""
    surfacePoints(uDegree::Int, vDegree::Int, uKnotVector, vKnotVector, controlPoints, uVector, vVector)

Compute B-spline surface: given the knotvectors and the degrees in 'u' and 'v' direction, the surface is evaluated at the evaluation points (uVector, vVector).

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
    uSpan = findSpan(nbasisFun, uVector, uKnotVector)
    Nu = basisFun(uSpan, uVector, uDegree, uKnotVector)

    # v-direction: determine the basis functions evaluated at vVector
    nbasisFun = length(vKnotVector) - vDegree - 1
    vSpan = findSpan(nbasisFun, vVector, vKnotVector)
    Nv = basisFun(vSpan, vVector, vDegree, vKnotVector)

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
