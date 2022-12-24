
"""
    surfacePoints(Patch::NURBSsurface, uEvalpoints, vEvalpoints)

Convenience function to compute points on a NURBSsurface.
"""
surfacePoints(Patch::NURBSsurface, uEvalpoints, vEvalpoints) = surfacePoints(
    Patch.uBasis.degree,
    Patch.vBasis.degree,
    Patch.uBasis.knotVec,
    Patch.vBasis.knotVec,
    Patch.controlPoints,
    uEvalpoints,
    vEvalpoints,
    Patch.weights,
)


"""
    surfacePoints(uDegree::Int, vDegree::Int, uKnotVector, vKnotVector, controlPoints, uVector, vVector, weights)

Compute NURBS surface: given the knotvectors and the degrees in 'u' and 'v' direction, the surface is evaluated at the evaluation points (uVector, vVector).

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

Note: the efficient evaluation via the B-spline basis is employed (no use of the naive evaluation of the NURBS basis).
"""
function surfacePoints(uDegree::Int, vDegree::Int, uKnotVector, vKnotVector, controlPoints, uVector, vVector, weights)

    # u-direction: determine the basis functions evaluated at uVector 
    nbasisFun = length(uKnotVector) - uDegree - 1
    uSpan = findSpan(nbasisFun, uVector, uKnotVector)
    Nu = basisFun(uSpan, uVector, uDegree, uKnotVector)

    # v-direction: determine the basis functions evaluated at vVector
    nbasisFun = length(vKnotVector) - vDegree - 1
    vSpan = findSpan(nbasisFun, vVector, vKnotVector)
    Nv = basisFun(vSpan, vVector, vDegree, vKnotVector)

    # intialize
    surface   = [SVector(0.0, 0.0, 0.0) for i in eachindex(uVector), j in eachindex(vVector)]
    normalize = [0.0 for i in eachindex(uVector), j in eachindex(vVector)]

    # determine the surface values
    for uPointInd in eachindex(uVector)

        uind = uSpan[uPointInd] - uDegree - 1

        for vPointInd in eachindex(vVector)

            for i in 1:(vDegree + 1)

                temp = SVector(0.0, 0.0, 0.0)
                normTemp = 0.0

                vind = vSpan[vPointInd] - vDegree + i - 1
                for k in 1:(uDegree + 1)

                    aux = Nu[uPointInd, k] * weights[uind + k, vind]

                    temp     += aux * controlPoints[uind + k, vind]
                    normTemp += aux
                end

                surface[uPointInd, vPointInd]   += Nv[vPointInd, i] * temp
                normalize[uPointInd, vPointInd] += Nv[vPointInd, i] * normTemp
            end

            surface[uPointInd, vPointInd] /= normalize[uPointInd, vPointInd]
        end
    end

    return surface
end
