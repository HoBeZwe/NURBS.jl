
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

    # Promote input types for initialization
    T = promote_type(eltype(eltype(controlPoints)), eltype(uVector))

    # determine the basis functions evaluated at uVector
    prealloc = preAlloc(degree(basis), uVector)
    basis(uVector, prealloc) # modifies prealloc

    # determine the curve values
    curve = fill(SVector{3,T}(0.0, 0.0, 0.0), length(uVector)) # initialize

    for (j, span) in enumerate(prealloc.spanVec)
        for ind in 1:(degree(basis) + 1)

            curve[j] += prealloc.B[j, ind] * controlPoints[span - degree(basis) + ind - 1]
        end
    end

    return curve
end


"""
    (curve::BsplineCurve)(uVector, k::Int)

Convenience function to compute points on all k derivatives of a B-spline curve.
"""
(curve::BsplineCurve)(uVector, k::Int) = curveDerivativesPoints(curve.basis, curve.controlPoints, uVector, k)


"""
    curvePoints(nbasisFun::Int, degree::Int, knotVector, controlPoints, uVector)

Compute a points on the k-th derivatives of a B-spline curve: given the 'knotVector', the 'controlPoints', and the 'degree', the curve is evaluated at the points given in 'uVector'.

Example for the controlPoints:

P1 = SVector(0.0, 0.0, 0.0)
P2 = SVector(0.1, 0.25, 0.0)
P3 = SVector(0.25, 0.3, 0.0)

controlPoints = [P1, P2, P3]
"""
function curveDerivativesPoints(basis::Basis, controlPoints, uVector, k::Int)

    # Promote input types for initialization
    T = promote_type(eltype(eltype(controlPoints)), eltype(uVector))

    # determine the basis functions evaluated at uVector
    prealloc = preAllocDer(degree(basis), uVector, k)
    basis(uVector, k, prealloc) # modifies prealloc

    # initialize
    curves = Vector{Vector{SVector{3,T}}}(undef, k + 1)
    for i in 1:(k + 1)
        curves[i] = fill(SVector{3,T}(0.0, 0.0, 0.0), length(uVector))
    end

    # determine the curve values
    for q in 1:(k + 1)
        for (j, span) in enumerate(prealloc.spanVec)
            for ind in 1:(degree(basis) + 1)

                curves[q][j] += prealloc.dersv[j, q, ind] * controlPoints[span - degree(basis) + ind - 1]
            end
        end
    end

    return curves
end
