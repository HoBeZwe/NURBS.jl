
"""
    (Patch::NURBSsurface)(uEvalpoints, vEvalpoints)

Convenience function to compute points on a NURBSsurface.
"""
(Patch::NURBSsurface)(uEvalpoints, vEvalpoints) = surfacePoints(
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


"""
    (Patch::NURBSsurface)(uEvalpoints, vEvalpoints, k::Int)

Convenience function to compute points on all k derivatives of a NURBSsurface.
"""
(Patch::NURBSsurface)(uEvalpoints, vEvalpoints, k::Int) = surfaceDerivativesPoints(
    Patch.uBasis.degree,
    Patch.vBasis.degree,
    Patch.uBasis.knotVec,
    Patch.vBasis.knotVec,
    Patch.controlPoints,
    uEvalpoints,
    vEvalpoints,
    Patch.weights,
    k,
)


"""
    surfaceDerivativesPoints(uDegree::Int, vDegree::Int, uKnotVector, vKnotVector, controlPoints, uVector, vVector, weights, k::Int)

Compute NURBS surface: given the knotvectors and the degrees in 'u' and 'v' direction, the surface is evaluated at the evaluation points (uVector, vVector).

Using (4.20) on page 136 of 'The NURBS Book'.

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

Returns a (k x k) matrix where each entry is a matrix of size (uKnotVector x vKnotVector): surfaces[q, p] is the matrix for the (q-1)-th derivative in u-direction and the (p-1)-th derivative in v-direction.

Note: the efficient evaluation via the B-spline basis is employed (no use of the naive evaluation of the NURBS basis).
"""
function surfaceDerivativesPoints(
    uDegree::Int, vDegree::Int, uKnotVector, vKnotVector, controlPoints, uVector, vVector, weights, k::Int
)

    # u-direction: determine the basis functions evaluated at uVector 
    nbasisFun = length(uKnotVector) - uDegree - 1
    uSpan = findSpan(nbasisFun, uVector, uKnotVector)
    Nu = derBasisFun(uSpan, uDegree, uVector, uKnotVector, k)

    # v-direction: determine the basis functions evaluated at vVector
    nbasisFun = length(vKnotVector) - vDegree - 1
    vSpan = findSpan(nbasisFun, vVector, vKnotVector)
    Nv = derBasisFun(vSpan, vDegree, vVector, vKnotVector, k)

    # intialize
    surfaces = [[SVector(0.0, 0.0, 0.0) for i in eachindex(uVector), j in eachindex(vVector)] for q in 1:(k + 1), p in (1:(k + 1))]
    w = [[0.0 for i in eachindex(uVector), j in eachindex(vVector)] for q in 1:(k + 1), p in 1:(k + 1)]


    # --- q-th derivative in u and p-th derivative in v
    for q in 0:k
        for p in 0:k

            # loop over u-points
            for uPointInd in eachindex(uVector)

                uind = uSpan[uPointInd] - uDegree - 1

                # loop over v-points
                for vPointInd in eachindex(vVector)

                    # derivatives of A and w
                    for ind in 1:(vDegree + 1)

                        temp = SVector(0.0, 0.0, 0.0)
                        normTemp = 0.0

                        vind = vSpan[vPointInd] - vDegree + ind - 1
                        for kind in 1:(uDegree + 1)

                            aux = Nu[uPointInd, q + 1, kind] * weights[uind + kind, vind]

                            temp     += aux * controlPoints[uind + kind, vind]
                            normTemp += aux
                        end

                        surfaces[q + 1, p + 1][uPointInd, vPointInd] += Nv[vPointInd, p + 1, ind] * temp
                        w[q + 1, p + 1][uPointInd, vPointInd]        += Nv[vPointInd, p + 1, ind] * normTemp
                    end

                    # second term
                    for i in 1:q
                        surfaces[q + 1, p + 1][uPointInd, vPointInd] -=
                            binomial(q, i) * w[i + 1, 1][uPointInd, vPointInd] * surfaces[q - i + 1, p + 1][uPointInd, vPointInd]
                    end

                    # third term
                    for j in 1:p
                        surfaces[q + 1, p + 1][uPointInd, vPointInd] -=
                            binomial(p, j) * w[1, j + 1][uPointInd, vPointInd] * surfaces[q + 1, p - j + 1][uPointInd, vPointInd]
                    end

                    # forth term
                    for i in 1:q
                        for j in 1:p
                            surfaces[q + 1, p + 1][uPointInd, vPointInd] -=
                                binomial(q, i) *
                                binomial(p, j) *
                                w[i + 1, j + 1][uPointInd, vPointInd] *
                                surfaces[q - i + 1, p - j + 1][uPointInd, vPointInd]
                        end
                    end

                    surfaces[q + 1, p + 1][uPointInd, vPointInd] /= w[1, 1][uPointInd, vPointInd]
                end
            end


        end
    end

    return surfaces
end
