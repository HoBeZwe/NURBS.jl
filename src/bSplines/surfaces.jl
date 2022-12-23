
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
