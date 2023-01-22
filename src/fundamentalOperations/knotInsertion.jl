
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
    oldMultIndices, oldMult = extendKnotVector!(knotVecOrig, degree, newParametricPoint, multiplicity)

    # --- modify control points
    extendControlPoints!(
        controlPointsOrig, knotVecCopy, degree, oldMultIndices.stop, newParametricPoint, multiplicity, oldMult, weights
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
    multiplicity + oldMult ≤ degree || error("The resulting multiplicity is larger than the polynomial degree.")

    # --- generate new knot vector
    toInsert = repeat([newParametricPoint], multiplicity + oldMult) # repeat the parametric point to the final multiplicity (old + new)
    splice!(knotVecOrig, oldMultIndices, toInsert)                  # insert the total parametric point in final multiplicity

    return oldMultIndices, oldMult
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

        for r in 1:multiplicity
            for i in 1:(degree - oldMult - r + 1)

                L = pos - degree + i + r - 1

                αᵢ = (uNew - knotVecOrig[L]) / (knotVecOrig[pos + i] - knotVecOrig[L])
                auxC[i] = αᵢ * auxW[i + 1] * auxC[i + 1] + (1.0 - αᵢ) * auxW[i] * auxC[i]
                auxW[i] = αᵢ * auxW[i + 1] + (1.0 - αᵢ) * auxW[i]
            end

            newControlPoints[r] = auxC[1] / auxW[1]
            newControlPoints[end - r + 1] = auxC[end - r] / auxW[end - r]

            newWeights[r] = auxW[1]
            newWeights[end - r + 1] = auxW[end - r]
        end

        # load remaining controlpoints
        for i in (multiplicity + oldMult + 1):(degree - oldMult - 1)
            newControlPoints[i] = auxC[i] / auxW[i]
            newWeights[i] = auxW[i]
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
    for i in (multiplicity + oldMult + 1):(degree - oldMult - 1)
        newControlPoints[i] = aux[i]
    end

    # splice new computed points into array (replacing the correct amount)
    splice!(controlPoints, (pos - degree + 1):(pos - oldMult - 1), newControlPoints)

    return nothing
end
