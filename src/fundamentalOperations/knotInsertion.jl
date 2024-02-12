
"""
    refine(C::Curve, newParametricPoints::Vector)

Convenience function to refine a curve.
"""
function refine(C::Curve, newParametricPoints::Vector)

    knotVecNew, controlPointsNew, weightsNew = refine(
        C.basis.knotVec, C.controlPoints, C.basis.degree, newParametricPoints, weights(C.basis)
    )

    return similarCurve(C, C.basis.degree, knotVecNew, controlPointsNew, weightsNew)
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
        insertKnot!(knotVecNew, controlPointsNew, degree, newPoint, 1, weightsNew)
    end

    return knotVecNew, controlPointsNew, weightsNew
end


"""
    insertKnot(C::Curve, newParametricPoint::Real, multiplicity::Int)

Convenience function to insert a knot into a curve.
"""
function insertKnot(C::Curve, newParametricPoint::Real, multiplicity::Int=1)

    knotVecNew, controlPointsNew, weightsNew = insertKnot(
        C.basis.knotVec, C.controlPoints, C.basis.degree, newParametricPoint, multiplicity, weights(C.basis)
    )

    return similarCurve(C, C.basis.degree, knotVecNew, controlPointsNew, weightsNew)
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
"""
function extendKnotVector!(knotVecOrig, degree::Int, newParametricPoint::Real, multiplicity::Int)

    # --- determine old multiplicity
    oldMultIndices = searchsorted(knotVecOrig, newParametricPoint)  # find existing entries that have the same value as the new parametric point
    oldMult = length(oldMultIndices)                                # number of existing entries with same value

    # --- check whether resulting mulitplicity makes sense
    newMult = multiplicity + oldMult
    if newMult > degree
        @info "The multiplicity of the inserted knot is limited to $(degree)."
        newMult = degree
    end

    # --- generate new knot vector
    toInsert = repeat([newParametricPoint], newMult) # repeat the parametric point to the final multiplicity (old + new)
    splice!(knotVecOrig, oldMultIndices, toInsert)                  # insert the total parametric point in final multiplicity

    return oldMultIndices, oldMult, newMult - oldMult
end


"""
    extendControlPoints!(controlPoints, knotVecOrig, degree::Int, pos::Int, uNew::Real, multiplicity::Int, oldMult::Int, weights)

Insert the new control points (and optionally the weights) corresponding to the new values in the knot vector.

Adaption of Algorithm A5.1 from 'The NURBS Book' p. 151.
"""
function extendControlPoints!(controlPoints, knotVecOrig, degree::Int, pos::Int, uNew::Real, multiplicity::Int, oldMult::Int, weights)

    T = eltype(controlPoints[1])

    # --- computation including weights
    if !isempty(weights)

        newControlPoints = [SVector{3,T}(0.0, 0.0, 0.0) for i in 1:(degree - oldMult + multiplicity - 1)]
        newWeights = [T(0.0) for i in 1:(degree - oldMult + multiplicity - 1)]

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
    newControlPoints = [SVector{3,T}(0.0, 0.0, 0.0) for i in 1:(degree - oldMult + multiplicity - 1)]

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
