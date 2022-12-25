
"""
    surfacePoints(Patch::BsplineSurface, uEvalpoints, vEvalpoints)

Convenience function to compute points on a B-spline surface.
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
                    temp += Nu[uPointInd, k] * controlPoints[uind + k, vind]
                end

                surface[uPointInd, vPointInd] += Nv[vPointInd, i] * temp
            end
        end
    end

    return surface
end


"""
    surfaceDerivativesPoints(Patch::BsplineSurface, uEvalpoints, vEvalpoints, k::Int)

Convenience function to compute points of the k-th derivatives of a B-spline surface.
"""
surfaceDerivativesPoints(Patch::BsplineSurface, uEvalpoints, vEvalpoints, k::Int) = surfaceDerivativesPoints(
    Patch.uBasis.degree,
    Patch.vBasis.degree,
    Patch.uBasis.knotVec,
    Patch.vBasis.knotVec,
    Patch.controlPoints,
    uEvalpoints,
    vEvalpoints,
    k,
)


"""
    surfaceDerivativesPoints(uDegree::Int, vDegree::Int, uKnotVector, vKnotVector, controlPoints, uVector, vVector, k::Int)

Compute B-spline surface and its derivatives: given the knotvectors and the degrees in 'u' and 'v' direction, the surface and its derivatives are evaluated at the evaluation points (uVector, vVector).

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

Note: not the most efficient implementation.
TODO: implement algorithm A.37 and A.38 of 'The Nurbs book'
"""
function surfaceDerivativesPoints(uDegree::Int, vDegree::Int, uKnotVector, vKnotVector, controlPoints, uVector, vVector, k::Int)

    # u-direction: determine the basis functions evaluated at uVector 
    nbasisFun = length(uKnotVector) - uDegree - 1
    uSpan = findSpan(nbasisFun, uVector, uKnotVector)
    Nu = derBasisFun(uSpan, uDegree, uVector, uKnotVector, k)

    # v-direction: determine the basis functions evaluated at vVector
    nbasisFun = length(vKnotVector) - vDegree - 1
    vSpan = findSpan(nbasisFun, vVector, vKnotVector)
    Nv = derBasisFun(vSpan, vDegree, vVector, vKnotVector, k)

    # intialize
    surfaces = [[SVector(0.0, 0.0, 0.0) for i in eachindex(uVector), j in eachindex(vVector)] for q in 1:(k + 1), p in 1:(k + 1)]

    # determine the surface derivative values
    for q in 1:(k + 1)
        for p in 1:(k + 1)
            surfaces[p, q] = surfaceDerivativesPointsUV(uDegree, vDegree, controlPoints, uVector, vVector, Nu, Nv, q, p, uSpan, vSpan) # derivatives along u and v
        end
    end

    return surfaces
end


"""
    surfaceDerivativesPointsUV(uDegree::Int, vDegree::Int, controlPoints, uVector, vVector, Nu, Nv, q::Int, p::Int, uSpan, vSpan)

Compute the q-th derivative along 'u' and the p-th derivative along 'v'.
"""
function surfaceDerivativesPointsUV(uDegree::Int, vDegree::Int, controlPoints, uVector, vVector, Nu, Nv, q::Int, p::Int, uSpan, vSpan)

    surfaces = [SVector(0.0, 0.0, 0.0) for i in eachindex(uVector), j in eachindex(vVector)]

    for uPointInd in eachindex(uVector)

        uind = uSpan[uPointInd] - uDegree - 1

        for vPointInd in eachindex(vVector)

            for i in 1:(vDegree + 1)

                temp = SVector(0.0, 0.0, 0.0)

                vind = vSpan[vPointInd] - vDegree + i - 1
                for kind in 1:(uDegree + 1)
                    temp += Nu[uPointInd, q, kind] * controlPoints[uind + kind, vind]
                end

                surfaces[uPointInd, vPointInd] += Nv[vPointInd, p, i] * temp
            end
        end
    end

    return surfaces
end
