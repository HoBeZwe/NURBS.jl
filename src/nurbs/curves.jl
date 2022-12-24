
"""
    curvePoints(curve::NURBScurve, uVector)

Convenience function to compute points on a NURBS curve.
"""
curvePoints(curve::NURBScurve, uVector) =
    curvePoints(curve.basis.degree, curve.basis.knotVec, curve.controlPoints, uVector, curve.basis.weights)


"""
    curvePoints(nbasisFun::Int, degree::Int, knotVector, controlPoints, uVector, weights)

Compute points on a NURBS curve: given the 'knotVector', the 'controlPoints', the 'degree', and the 'weights', the curve is evaluated at the points given in 'uVector'.

Example for the controlPoints:

P1 = SVector(0.0, 0.0, 0.0)
P2 = SVector(0.1, 0.25, 0.0)
P3 = SVector(0.25, 0.3, 0.0)

controlPoints = [P1, P2, P3]

Note: the efficient evaluation via the B-spline basis is employed (no use of the naive evaluation of the NURBS basis).
"""
function curvePoints(degree::Int, knotVector, controlPoints, uVector, weights)

    # the number of basis functions is determined by the number of knot vectors and the degree
    nbasisFun = length(knotVector) - degree - 1

    # determine the basis functions evaluated at uVector
    spans = findSpan(nbasisFun, uVector, knotVector)
    N = basisFun(spans, uVector, degree, knotVector)

    # determine the curve values
    curve     = [SVector(0.0, 0.0, 0.0) for i in eachindex(uVector)] # initialize
    normalize = [0.0 for i in eachindex(uVector)]

    for (j, span) in enumerate(spans)
        for ind in 1:(degree + 1)

            index = span - degree + ind - 1

            aux = N[j, ind] * weights[index]

            curve[j]     += aux * controlPoints[index]
            normalize[j] += aux
        end

        curve[j] /= normalize[j]

        normalize[j] == 0.0 && error("division by zero!") # TODO: properly handle this
    end

    return curve
end
