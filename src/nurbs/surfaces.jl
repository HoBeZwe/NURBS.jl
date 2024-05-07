
"""
    (Patch::NURBSsurface)(uEvalpoints, vEvalpoints)

Convenience function to compute points on a NURBSsurface.
"""
(Patch::NURBSsurface)(uEvalpoints, vEvalpoints) =
    surfacePoints(Patch.uBasis, Patch.vBasis, Patch.controlPoints, uEvalpoints, vEvalpoints, Patch.weights)


"""
    surfacePoints(uBasis::Basis, vBasis::Basis, controlPoints, uVector, vVector, weights)

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
function surfacePoints(uBasis::Basis, vBasis::Basis, controlPoints, uVector, vVector, weights)

    # Promote input types for initialization
    T = promote_type(eltype.(controlPoints[:])..., eltype(uVector), eltype(vVector), eltype(weights))

    # u-direction: determine the basis functions evaluated at uVector
    nbasisFun = numBasisFunctions(uBasis)
    uSpan = findSpan(nbasisFun, uVector, uBasis.knotVec, uBasis.degree)
    Nu = basisFun(uSpan, uVector, uBasis)

    # v-direction: determine the basis functions evaluated at vVector
    nbasisFun = numBasisFunctions(vBasis)
    vSpan = findSpan(nbasisFun, vVector, vBasis.knotVec, vBasis.degree)
    Nv = basisFun(vSpan, vVector, vBasis)

    # intialize
    surface   = [SVector{3,T}(0.0, 0.0, 0.0) for i in eachindex(uVector), j in eachindex(vVector)]
    normalize = [T(0.0) for i in eachindex(uVector), j in eachindex(vVector)]

    # determine the surface values
    for uPointInd in eachindex(uVector)

        uind = uSpan[uPointInd] - uBasis.degree - 1

        for vPointInd in eachindex(vVector)

            for i in 1:(vBasis.degree + 1)

                temp = SVector{3,T}(0.0, 0.0, 0.0)
                normTemp = T(0.0)

                vind = vSpan[vPointInd] - vBasis.degree + i - 1
                for k in 1:(uBasis.degree + 1)

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


struct pAllocNURBSsuface{T<:Real,F<:Int,L}
    preallocU::pAllocDer{T,F,L}
    preallocV::pAllocDer{T,F,L}
    surfaces::Matrix{Matrix{SVector{3,T}}}
    w::Matrix{Matrix{T}}
end


"""
    (Patch::NURBSsurface)(uEvalpoints, vEvalpoints, k::Int, prealloc)

Convenience function to compute points on all k derivatives of a NURBSsurface, for preallocated memory.
"""
(Patch::NURBSsurface)(uEvalpoints, vEvalpoints, k::Int, prealloc::pAllocNURBSsuface) = surfaceDerivativesPoints!(
    prealloc,
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

Allocate memory and call surfaceDerivativesPoints!
"""
function surfaceDerivativesPoints(
    uDegree::Int, vDegree::Int, uKnotVector, vKnotVector, controlPoints, uVector, vVector, weights, k::Int
)
    prealloc = preAllocNURBSsurface(uDegree, vDegree, uVector, vVector, k)

    return surfaceDerivativesPoints!(prealloc, uDegree, vDegree, uKnotVector, vKnotVector, controlPoints, uVector, vVector, weights, k)
end


"""
    preAllocNURBSsurface(uDegree::Int, vDegree::Int, uVector, vVector, k::Int)

Allocate all memory for surfaceDerivativesPoints!
"""
function preAllocNURBSsurface(uDegree::Int, vDegree::Int, uVector, vVector, k::Int)

    # Promote input types for initialization
    T = promote_type(eltype(uVector), eltype(vVector))

    preallocU = preAllocDer(uDegree, uVector, k)
    preallocV = preAllocDer(vDegree, vVector, k)

    surfaces = [
        [SVector{3,T}(0.0, 0.0, 0.0) for i in eachindex(uVector), j in eachindex(vVector)] for q in 1:(k + 1), p in (1:(k + 1))
    ]
    w = [[T(0.0) for i in eachindex(uVector), j in eachindex(vVector)] for q in 1:(k + 1), p in 1:(k + 1)]

    return pAllocNURBSsuface(preallocU, preallocV, surfaces, w)
end



"""
    surfaceDerivativesPoints!(prealloc, uDegree::Int, vDegree::Int, uKnotVector, vKnotVector, controlPoints, uVector, vVector, weights, k::Int)

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
function surfaceDerivativesPoints!(
    prealloc::pAllocNURBSsuface, uDegree::Int, vDegree::Int, uKnotVector, vKnotVector, controlPoints, uVector, vVector, weights, k::Int
)

    # Promote input types for initialization
    T = promote_type(
        eltype(uKnotVector), eltype(vKnotVector), eltype.(controlPoints[:])..., eltype(uVector), eltype(vVector), eltype(weights)
    )

    preallocU = prealloc.preallocU
    preallocV = prealloc.preallocV
    surfaces = prealloc.surfaces
    w = prealloc.w

    # u-direction: determine the basis functions evaluated at uVector
    nbasisFun = length(uKnotVector) - uDegree - 1
    uSpan = findSpan!(preallocU.spanVec, nbasisFun, uVector, uKnotVector, uDegree)
    Nu = derBasisFun!(preallocU, uSpan, uDegree, uVector, uKnotVector, k)

    # v-direction: determine the basis functions evaluated at vVector
    nbasisFun = length(vKnotVector) - vDegree - 1
    vSpan = findSpan!(preallocV.spanVec, nbasisFun, vVector, vKnotVector, vDegree)
    Nv = derBasisFun!(preallocV, vSpan, vDegree, vVector, vKnotVector, k)

    # initialize
    for i in eachindex(w)
        for j in eachindex(w[1])
            w[i][j] = 0.0
        end
    end

    for i in eachindex(surfaces)
        for j in eachindex(surfaces[1])
            surfaces[i][j] = SVector{3,T}(0.0, 0.0, 0.0)
        end
    end

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

                        temp = SVector{3,T}(0.0, 0.0, 0.0)
                        normTemp = T(0.0)

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
