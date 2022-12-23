
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
