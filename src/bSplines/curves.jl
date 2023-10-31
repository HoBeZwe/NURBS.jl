
"""
    (curve::BsplineCurve)(uVector)

Convenience function to compute points on a B-spline curve.
"""
(curve::BsplineCurve)(uVector) = curvePoints(curve.basis, curve.controlPoints, uVector)


"""
    curvePoints(basis::Basis, controlPoints, uVector)

Compute a 1D B-spline curve: given the 'knotVector', the 'controlPoints', and the 'degree', the curve is evaluated at the points given in 'uVector'.

Example for the controlPoints:
    
P1 = SVector(0.0, 0.0, 0.0)
P2 = SVector(0.1, 0.25, 0.0)
P3 = SVector(0.25, 0.3, 0.0)

controlPoints = [P1, P2, P3]
"""
function curvePoints(basis::Basis, controlPoints, uVector)

    # the number of basis functions is determined by the number of knot vectors and the degree
    nbasisFun = numBasisFunctions(basis)

    # determine the basis functions evaluated at uVector
    spans = findSpan(nbasisFun, uVector, basis.knotVec, basis.degree)
    N = basisFun(spans, uVector, basis)

    # determine the curve values
    curve = [SVector(0.0, 0.0, 0.0) for i in eachindex(uVector)] # initialize

    for (j, span) in enumerate(spans)
        for ind in 1:(basis.degree + 1)

            curve[j] += N[j, ind] * controlPoints[span - basis.degree + ind - 1]
        end
    end

    return curve
end


"""
    (curve::BsplineCurve)(uVector, k::Int)

Convenience function to compute points on all k derivatives of a B-spline curve.
"""
(curve::BsplineCurve)(uVector, k::Int) =
    curveDerivativesPoints(curve.basis.degree, curve.basis.knotVec, curve.controlPoints, uVector, k)


"""
    curvePoints(nbasisFun::Int, degree::Int, knotVector, controlPoints, uVector)

Compute a points on the k-th derivatives of a B-spline curve: given the 'knotVector', the 'controlPoints', and the 'degree', the curve is evaluated at the points given in 'uVector'.

Example for the controlPoints:
    
P1 = SVector(0.0, 0.0, 0.0)
P2 = SVector(0.1, 0.25, 0.0)
P3 = SVector(0.25, 0.3, 0.0)

controlPoints = [P1, P2, P3]
"""
function curveDerivativesPoints(degree::Int, knotVector, controlPoints, uVector, k::Int)

    # the number of basis functions is determined by the number of knot vectors and the degree
    nbasisFun = length(knotVector) - degree - 1

    # determine the basis functions evaluated at uVector
    spans = findSpan(nbasisFun, uVector, knotVector, degree)
    N = derBasisFun(spans, degree, uVector, knotVector, k)

    # determine the curve values
    #curve = [SVector(0.0, 0.0, 0.0) for i in eachindex(uVector)] # initialize
    curves = [[SVector(0.0, 0.0, 0.0) for i in eachindex(uVector)] for j in 1:(k + 1)]

    for q in 1:(k + 1)
        for (j, span) in enumerate(spans)
            for ind in 1:(degree + 1)

                curves[q][j] += N[j, q, ind] * controlPoints[span - degree + ind - 1]
            end
        end
    end

    return curves
end
