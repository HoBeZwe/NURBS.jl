
"""
    (Patch::BsplineSurface)(uEvalpoints, vEvalpoints)

Convenience function to compute points on a B-spline surface.
"""
(Patch::BsplineSurface)(uEvalpoints, vEvalpoints) =
    surfacePoints(Patch.uBasis, Patch.vBasis, Patch.controlPoints, uEvalpoints, vEvalpoints)


"""
    surfacePoints(uBasis::Basis, vBasis::Basis, controlPoints, uVector, vVector)

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
function surfacePoints(uBasis::Basis, vBasis::Basis, controlPoints, uVector, vVector)

    # Promote input types for initialization
    T = promote_type(eltype(eltype(controlPoints)), eltype(uVector), eltype(vVector))

    # u-direction: determine the basis functions evaluated at uVector
    preallocU = preAlloc(degree(uBasis), uVector)
    uBasis(uVector, preallocU) # prealloc is modified

    # v-direction: determine the basis functions evaluated at vVector
    preallocV = preAlloc(degree(vBasis), vVector)
    vBasis(vVector, preallocV) # prealloc is modified

    # intialize
    surface = fill(SVector{3,T}(0.0, 0.0, 0.0), length(uVector), length(vVector))

    # determine the surface values
    for uPointInd in eachindex(uVector)

        uind = preallocU.spanVec[uPointInd] - degree(uBasis) - 1

        for vPointInd in eachindex(vVector)

            for i in 1:(degree(vBasis) + 1)

                temp = SVector{3,T}(0.0, 0.0, 0.0)

                vind = preallocV.spanVec[vPointInd] - degree(vBasis) + i - 1
                for k in 1:(degree(uBasis) + 1)
                    temp += preallocU.B[uPointInd, k] * controlPoints[uind + k, vind]
                end

                surface[uPointInd, vPointInd] += preallocV.B[vPointInd, i] * temp
            end
        end
    end

    return surface
end


"""
    (Patch::BsplineSurface)(uEvalpoints, vEvalpoints, k::Int)

Convenience function to compute points on all k derivatives of a B-spline surface.
"""
(Patch::BsplineSurface)(uEvalpoints, vEvalpoints, k::Int) =
    surfaceDerivativesPoints(Patch.uBasis, Patch.vBasis, Patch.controlPoints, uEvalpoints, vEvalpoints, k)


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

Returns a (k x k) matrix where each entry is a matrix of size (uKnotVector x vKnotVector): surfaces[q, p] is the matrix for the (q-1)-th derivative in u-direction and the (p-1)-th derivative in v-direction.

Note: not the most efficient implementation.
TODO: implement algorithm A.37 and A.38 of 'The Nurbs book'
"""
function surfaceDerivativesPoints(uBasis::Basis, vBasis::Basis, controlPoints, uVector, vVector, k::Int)

    # Promote input types for initialization
    T = promote_type(eltype(eltype(controlPoints)), eltype(uVector), eltype(vVector))

    # u-direction: determine the basis functions evaluated at uVector
    preallocU = preAllocDer(degree(uBasis), uVector, k)
    uBasis(uVector, k, preallocU) # prealloc is modified

    # v-direction: determine the basis functions evaluated at vVector
    preallocV = preAllocDer(degree(vBasis), vVector, k)
    vBasis(vVector, k, preallocV) # prealloc is modified

    # intialize
    kp = k + 1
    lu = length(uVector)
    lv = length(vVector)
    surfaces = Matrix{Matrix{SVector{3,T}}}(undef, kp, kp)
    for i in 1:kp, j in 1:kp
        surfaces[i, j] = fill(SVector{3,T}(0.0, 0.0, 0.0), lu, lv)
    end

    # determine the surface derivative values
    for q in 1:kp
        for p in 1:kp
            surfaces[q, p] = surfaceDerivativesPointsUV(
                degree(uBasis),
                degree(vBasis),
                controlPoints,
                uVector,
                vVector,
                preallocU.dersv,
                preallocV.dersv,
                q,
                p,
                preallocU.spanVec,
                preallocV.spanVec,
            ) # derivatives along u and v
        end
    end

    return surfaces
end


"""
    surfaceDerivativesPointsUV(uDegree::Int, vDegree::Int, controlPoints, uVector, vVector, Nu, Nv, q::Int, p::Int, uSpan, vSpan)

Compute the q-th derivative along 'u' and the p-th derivative along 'v'.
"""
function surfaceDerivativesPointsUV(uDegree::Int, vDegree::Int, controlPoints, uVector, vVector, Nu, Nv, q::Int, p::Int, uSpan, vSpan)

    # Promote input types for initialization
    T = promote_type(eltype(eltype(controlPoints)), eltype(uVector), eltype(vVector))

    surfaces = fill(SVector{3,T}(0.0, 0.0, 0.0), length(uVector), length(uVector))

    for uPointInd in eachindex(uVector)

        uind = uSpan[uPointInd] - uDegree - 1

        for vPointInd in eachindex(vVector)

            for i in 1:(vDegree + 1)

                temp = SVector{3,T}(0.0, 0.0, 0.0)

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
