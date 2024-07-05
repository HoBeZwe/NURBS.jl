
"""
    checkRange(point<:Real)

Check whether points is greater 0 and smaller 1.
"""
function checkRange(point::Real)

    (point > 0.0 && point < 1.0) || error("A parametric point is outside ]0,1[.")

    return nothing
end


"""
    refine(C::Curve, newParametricPoints::Vector)

Convenience function to refine a curve.
"""
function refine(C::Curve, newParametricPoints::Vector)

    knotVecNew, controlPointsNew, weightsNew = refine(
        C.basis.knotVec, C.controlPoints, degree(C), newParametricPoints, weights(C.basis)
    )

    return similarCurve(C, degree(C), knotVecNew, controlPointsNew, weightsNew)
end


"""
    refine(kVec::Vector, controlPts::Vector, degree::Int, newParametricPoints::Vector, weights=[])

Refine a curve by inserting new parametric points into the curve's knot vector.

NOTE: There are more efficient ways to do this. See, e.g., 'The NURBS Book' p. 162.
"""
function refine(kVec::Vector, controlPts::Vector, degree::Int, newParametricPoints::Vector, weights=[])

    # --- make copies
    knotVecNew       = deepcopy(kVec)
    controlPointsNew = deepcopy(controlPts)
    weightsNew       = deepcopy(weights)

    # --- modify copies
    for (i, newPoint) in enumerate(newParametricPoints)
        checkRange(newPoint)
        insertKnot!(knotVecNew, controlPointsNew, degree, newPoint, 1, weightsNew)
    end

    return knotVecNew, controlPointsNew, weightsNew
end


"""
    refine(S::Surface; U=[], V=[])

NOTE: There are more efficient ways to do this.
"""
function refine(S::Surface; U=[], V=[])

    for u in U
        checkRange(u)
        S = insertKnotU(S, u)
    end

    for v in V
        checkRange(v)
        S = insertKnotV(S, v)
    end

    return S
end


"""
    insertKnot(C::Curve, newParametricPoint::Real, multiplicity::Int)

Convenience function to insert a knot into a curve.
"""
function insertKnot(C::Curve, newParametricPoint::Real, multiplicity::Int=1)

    knotVecNew, controlPointsNew, weightsNew = insertKnot(
        C.basis.knotVec, C.controlPoints, degree(C), newParametricPoint, multiplicity, weights(C.basis)
    )

    return similarCurve(C, degree(C), knotVecNew, controlPointsNew, weightsNew)
end


"""
    insertKnotU(S::Surface, newParametricPoint::Real, multiplicity::Int=1)

Insert a knot in u-direction with the given multiplicity.

TODO: the computation of the alphas is redundant.
"""
function insertKnotU(S::Surface, newParametricPoint::Real, multiplicity::Int=1)

    degr = degree(S.uBasis)
    knotVecOrig = S.uBasis.knotVec

    # --- modify knot vector
    knotVecMod = deepcopy(knotVecOrig)
    oldMultIndices, oldMult, limitedMult = extendKnotVector!(knotVecMod, degr, newParametricPoint, multiplicity)

    # --- modify control points
    N = numBasisFunctions(knotVecMod, degr)
    M = size(S.controlPoints, 2)

    wMatNew = similar(weights(S), N, M)
    cPtsNew = similar(S.controlPoints, N, M)

    for i in eachindex(S.controlPoints[1, :])

        ctrlPtsMody = deepcopy(S.controlPoints[:, i])
        weightsMody = deepcopy(weights(S, :, i))

        extendControlPoints!(
            ctrlPtsMody, knotVecOrig, degr, oldMultIndices.stop, newParametricPoint, limitedMult, oldMult, weightsMody
        )

        isempty(weightsMody) || (wMatNew[:, i] = weightsMody)
        cPtsNew[:, i] = ctrlPtsMody
    end

    return similarSurface(S, degr, degree(S.vBasis), knotVecMod, S.vBasis.knotVec, cPtsNew, wMatNew)
end


"""
    insertKnotV(S::Surface, newParametricPoint::Real, multiplicity::Int=1)

Insert a knot in v-direction with the given multiplicity.

TODO: the computation of the alphas is redundant.
"""
function insertKnotV(S::Surface, newParametricPoint::Real, multiplicity::Int=1)

    degr = degree(S.vBasis)
    knotVecOrig = S.vBasis.knotVec

    # --- modify knot vector
    knotVecMod = deepcopy(knotVecOrig)
    oldMultIndices, oldMult, limitedMult = extendKnotVector!(knotVecMod, degr, newParametricPoint, multiplicity)

    # --- modify control points
    N = numBasisFunctions(knotVecMod, degr)
    M = size(S.controlPoints, 1)

    wMatNew = similar(weights(S), M, N)
    cPtsNew = similar(S.controlPoints, M, N)

    for i in eachindex(S.controlPoints[:, 1])

        ctrlPtsMody = deepcopy(S.controlPoints[i, :])
        weightsMody = deepcopy(weights(S, i, :))

        extendControlPoints!(
            ctrlPtsMody, knotVecOrig, degr, oldMultIndices.stop, newParametricPoint, limitedMult, oldMult, weightsMody
        )

        isempty(weightsMody) || (wMatNew[i, :] = weightsMody)
        cPtsNew[i, :] = ctrlPtsMody
    end

    return similarSurface(S, degree(S.uBasis), degr, S.uBasis.knotVec, knotVecMod, cPtsNew, wMatNew)
end


"""
    insertKnot(knotVecOrig, controlPointsOrig, p::Int, newParametricPoint::Real, multiplicity::Int)

Insert a new value with a given mulitplicity (insert the value multiple times) into a knot vector for the polynomial degree 'degree'. Return the resutling knot vector and the new control points.

Adaption of Algorithm A5.1 from 'The NURBS Book' p. 151.

The provided knot vector and control-point vector are NOT modified.
"""
function insertKnot(knotVecOrig, controlPointsOrig, degree::Int, newParametricPoint::Real, multiplicity::Int, weights=[])

    # --- make copies
    knotVecNew       = deepcopy(knotVecOrig)
    controlPointsNew = deepcopy(controlPointsOrig)
    weightsNew       = deepcopy(weights)

    # --- modify copies
    insertKnot!(knotVecNew, controlPointsNew, degree, newParametricPoint, multiplicity, weightsNew)

    return knotVecNew, controlPointsNew, weightsNew
end


"""
    insertKnot!(knotVecOrig, controlPointsOrig, p::Int, newParametricPoint::Real, multiplicity::Int)

Insert a new value with a given mulitplicity (insert the value multiple times) into a knot vector for the polynomial degree 'degree'. Return the resutling knot vector and the new control points.

Adaption of Algorithm A5.1 from 'The NURBS Book' p. 151.

The provided knot vector and control-point vector are modified.
"""
function insertKnot!(knotVecOrig, controlPointsOrig, degree::Int, newParametricPoint::Real, multiplicity::Int, weights=[])

    # --- is provide knot vector a valid one?
    isValidKnotVector!(knotVecOrig)

    # --- modify knot vector
    knotVecCopy = deepcopy(knotVecOrig)
    oldMultIndices, oldMult, limitedMult = extendKnotVector!(knotVecOrig, degree, newParametricPoint, multiplicity)

    # --- modify control points
    extendControlPoints!(
        controlPointsOrig, knotVecCopy, degree, oldMultIndices.stop, newParametricPoint, limitedMult, oldMult, weights
    )

    return nothing
end


"""
    extendKnotVector!(knotVecOrig, newParametricPoint::Real, multiplicity::Int)

Insert the new value with a given mulitplicity (insert the point multiple times).

Modifies knotVecOrig.
"""
function extendKnotVector!(knotVecOrig, degree::Int, newParametricPoint::Real, multiplicity::Int)

    # --- determine old multiplicity
    oldMultIndices = searchsorted(knotVecOrig, newParametricPoint)  # find existing entries that have the same value as the new parametric point
    oldMult = length(oldMultIndices)                                # number of existing entries with same value

    # --- check whether resulting mulitplicity makes sense
    newMult = multiplicity + oldMult
    newMult = limitMultiplicity(newMult, degree)

    # --- generate new knot vector
    toInsert = repeat([newParametricPoint], newMult) # repeat the parametric point to the final multiplicity (old + new)
    splice!(knotVecOrig, oldMultIndices, toInsert)   # insert the total parametric point in final multiplicity

    return oldMultIndices, oldMult, newMult - oldMult
end


"""
    limitMultiplicity(newMult::Int, degree::Int)

Limit the multiplicity to the maximum allowed value.
"""
function limitMultiplicity(newMult::Int, degree::Int)

    if degree == 0
        newMult > 1 && @info "The multiplicity of the inserted knot is limited to 1."
        return 1
    end

    if newMult > degree
        @info "The multiplicity of the inserted knot is limited to $(degree)."
        return degree
    end

    return newMult
end


"""
    extendControlPoints!(controlPoints, knotVecOrig, degree::Int, pos::Int, uNew::Real, multiplicity::Int, oldMult::Int, weights)

Insert the new control points (and optionally the weights) corresponding to the new values in the knot vector 'multiplicity' times.

Adaption of Algorithm A5.1 from 'The NURBS Book' p. 151.

Modifies controlPoints and weights.
"""
function extendControlPoints!(
    controlPoints, knotVecOrig, degree::Int, pos::Int, uNew::Real, multiplicity::Int, oldMult::Int, weights=[]
)

    degree == 0 && return extendControlPoints!(controlPoints, pos, multiplicity, oldMult)

    # --- computation including weights
    if !isempty(weights)

        newControlPoints = similar(controlPoints, degree - oldMult + multiplicity - 1) #[SVector{3,T}(0.0, 0.0, 0.0) for i in 1:(degree - oldMult + multiplicity - 1)]
        newWeights = similar(weights, degree - oldMult + multiplicity - 1) #[T(0.0) for i in 1:(degree - oldMult + multiplicity - 1)]

        auxC = controlPoints[(pos - degree):(pos - oldMult)]
        auxW = weights[(pos - degree):(pos - oldMult)]
        auxC .*= auxW

        for r in 1:multiplicity
            for i in 1:(degree - oldMult - r + 1)

                L = pos - degree + i + r - 1

                αᵢ = (uNew - knotVecOrig[L]) / (knotVecOrig[pos + i] - knotVecOrig[L])
                auxC[i] = αᵢ * auxC[i + 1] + (1.0 - αᵢ) * auxC[i]
                auxW[i] = αᵢ * auxW[i + 1] + (1.0 - αᵢ) * auxW[i]
            end

            newControlPoints[r] = auxC[1] / auxW[1]
            newControlPoints[end - r + 1] = auxC[end - r] / auxW[end - r]

            newWeights[r] = auxW[1]
            newWeights[end - r + 1] = auxW[end - r]
        end

        # load remaining controlpoints
        for i in 1:(degree - oldMult - 1)
            newControlPoints[i + multiplicity - 1] = auxC[i] / auxW[i]
            newWeights[i + multiplicity - 1] = auxW[i]
        end

        # splice new computed points into array (replacing the correct amount)
        splice!(controlPoints, (pos - degree + 1):(pos - oldMult - 1), newControlPoints)
        splice!(weights, (pos - degree + 1):(pos - oldMult - 1), newWeights)

        return nothing
    end

    # --- computation when there are no weights
    newControlPoints = similar(controlPoints, degree - oldMult + multiplicity - 1) #[SVector{3,T}(0.0, 0.0, 0.0) for i in 1:(degree - oldMult + multiplicity - 1)]

    aux = controlPoints[(pos - degree):(pos - oldMult)]
    for r in 1:multiplicity
        for i in 1:(degree - oldMult - r + 1)

            L = pos - degree + i + r - 1

            αᵢ = (uNew - knotVecOrig[L]) / (knotVecOrig[pos + i] - knotVecOrig[L])
            aux[i] = αᵢ * aux[i + 1] + (1.0 - αᵢ) * aux[i]
        end
        newControlPoints[r] = aux[1]
        newControlPoints[end - r + 1] = aux[end - r]
    end

    # load remaining controlpoints
    for i in 1:(degree - oldMult - 1)
        newControlPoints[i + multiplicity - 1] = aux[i]
    end

    # splice new computed points into array (replacing the correct amount)
    splice!(controlPoints, (pos - degree + 1):(pos - oldMult - 1), newControlPoints)

    return nothing
end

"""
Lowest order polynomial degree = 0.

Can multiplicity be maximum 1 and/or oldMult maximum 0? => simplify below algorithm further.
"""
function extendControlPoints!(controlPoints, pos::Int, multiplicity::Int, oldMult::Int)

    multiplicity == 0 && return nothing

    # splice new points into array 
    splice!(controlPoints, (pos + 1):pos, [controlPoints[pos]])

    return nothing
end
