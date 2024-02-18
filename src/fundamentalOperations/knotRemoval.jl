
"""
    removeKnot(C::Curve, pointToRemove::Real, multiplicity::Int)

Remove a knot from a curve ´multiplicity´ times.
"""
function removeKnot(C::Curve, pointToRemove::Real, multiplicity::Int)

    knotVecNew, controlPointsNew, weightsNew = removeKnot(
        C.basis.knotVec, C.controlPoints, C.basis.degree, pointToRemove, multiplicity, weights(C.basis)
    )

    return similarCurve(C, C.basis.degree, knotVecNew, controlPointsNew, weightsNew)
end


"""
    removeKnotU(S::Surface, pointToRemove::Real, multiplicity::Int)

Remove a knot in u-direction 'multiplicity' times.

TODO: replace by more efficient algorithm of Tiller.
"""
function removeKnotU(S::Surface, pointToRemove::Real, multiplicity::Int)

    degree = S.uBasis.degree
    knotVec = S.uBasis.knotVec

    # --- analyze the knot vector: how many times is the knot actually contained?
    oldMultIndices, oldMult, limitedMult = trimability(knotVec, pointToRemove, multiplicity)

    # --- remove knots from the knot vector
    knotVecMod = deepcopy(knotVec)
    deleteat!(knotVecMod, oldMultIndices[1:limitedMult])

    # --- modify control points
    N = numBasisFunctions(knotVecMod, degree)
    M = size(S.controlPoints, 2)

    wMatNew = similar(weights(S), N, M)
    cPtsNew = similar(S.controlPoints, N, M)

    for i in eachindex(S.controlPoints[1, :])

        ctrlPtsMody = deepcopy(S.controlPoints[:, i])
        weightsMody = deepcopy(weights(S, :, i))

        actualRemoved = trimControlPoints!(
            ctrlPtsMody, knotVec, degree, oldMultIndices.stop, pointToRemove, limitedMult, oldMult, weightsMody
        )

        if actualRemoved < limitedMult
            @warn "The knot $(pointToRemove) can not be removed $multiplicity times."
            return S
        end

        isempty(weightsMody) || (wMatNew[:, i] = weightsMody)
        cPtsNew[:, i] = ctrlPtsMody
    end

    return similarSurface(S, degree, S.vBasis.degree, knotVecMod, S.vBasis.knotVec, cPtsNew, wMatNew)
end


"""
    removeKnotV(S::Surface, pointToRemove::Real, multiplicity::Int)

Remove a knot in v-direction 'multiplicity' times.

TODO: replace by more efficient algorithm of Tiller.
"""
function removeKnotV(S::Surface, pointToRemove::Real, multiplicity::Int)

    degree = S.vBasis.degree
    knotVec = S.vBasis.knotVec

    # --- analyze the knot vector: how many times is the knot actually contained?
    oldMultIndices, oldMult, limitedMult = trimability(knotVec, pointToRemove, multiplicity)

    # --- remove knots from the knot vector
    knotVecMod = deepcopy(knotVec)
    deleteat!(knotVecMod, oldMultIndices[1:limitedMult])

    # --- modify control points
    N = numBasisFunctions(knotVecMod, degree)
    M = size(S.controlPoints, 1)

    wMatNew = similar(weights(S), M, N)
    cPtsNew = similar(S.controlPoints, M, N)

    for i in eachindex(S.controlPoints[:, 1])

        ctrlPtsMody = deepcopy(S.controlPoints[i, :])
        weightsMody = deepcopy(weights(S, i, :))

        actualRemoved = trimControlPoints!(
            ctrlPtsMody, knotVec, degree, oldMultIndices.stop, pointToRemove, limitedMult, oldMult, weightsMody
        )

        if actualRemoved < limitedMult
            @warn "The knot $(pointToRemove) can not be removed $multiplicity times."
            return S
        end

        isempty(weightsMody) || (wMatNew[i, :] = weightsMody)
        cPtsNew[i, :] = ctrlPtsMody
    end

    return similarSurface(S, S.uBasis.degree, degree, S.uBasis.knotVec, knotVecMod, cPtsNew, wMatNew)
end


"""
    removeKnot(knotVec, controlPoints, degree::Int, pointToRemove::Real, multiplicity::Int, weights)

Remove a knot 'multiplicity' times.
"""
function removeKnot(knotVec, controlPoints, degree::Int, pointToRemove::Real, multiplicity::Int, weights=[])

    # --- make copies
    knotVecMod = deepcopy(knotVec)
    ctrlPtsMod = deepcopy(controlPoints)
    weightsMod = deepcopy(weights)

    # --- modify copies
    removeKnot!(knotVecMod, ctrlPtsMod, degree, pointToRemove, multiplicity, weightsMod)

    return knotVecMod, ctrlPtsMod, weightsMod
end


"""
    removeKnot!(knotVec, controlPoints, degree::Int, pointToRemove::Real, multiplicity::Int, weights)

Modifies solely 'knotVec', 'controlPoints', and 'weights'.
"""
function removeKnot!(knotVec, controlPoints, degree::Int, pointToRemove::Real, multiplicity::Int, weights=[])

    # --- is provide knot vector a valid one?
    isValidKnotVector!(knotVec)

    # --- analyze the knot vector: how many times is the knot actually contained?
    oldMultIndices, oldMult, limitedMult = trimability(knotVec, pointToRemove, multiplicity)

    # --- modify control points
    actualRemoved = trimControlPoints!(
        controlPoints, knotVec, degree, oldMultIndices.stop, pointToRemove, limitedMult, oldMult, weights
    )

    # --- remove actually removable knots from the knot vector
    deleteat!(knotVec, oldMultIndices[1:actualRemoved])

    return nothing
end


"""
    trimability(knotVec, pointToRemove::Real, multiplicity::Int)

Analyze the knot vector: how many times is the knot actually contained, at which positions, and how many times.
"""
function trimability(knotVec, pointToRemove::Real, multiplicity::Int)

    # --- interior knot?
    if iszero(pointToRemove) || pointToRemove == 1.0
        error("The knot has to be an interior knot.") # assumption according to the algorithm p. 179       
    end

    oldMultIndices = searchsorted(knotVec, pointToRemove)  # find existing entries that have the same value as the point to be removed
    oldMult = length(oldMultIndices)                       # number of existing entries with same value

    # --- are the knots removable?
    limitedMult = multiplicity
    if multiplicity > oldMult
        if iszero(oldMult)
            @warn "The requested knot is not contained in the knot vector."
            return oldMultIndices, oldMult, 0
        end

        @warn "The knot is only contained $(oldMult) times."
        limitedMult = oldMult
    end

    return oldMultIndices, oldMult, limitedMult
end


"""
    trimControlPoints!(cPts, kVecOri, degree::Int, pos::Int, pointToRemove::Real, limitedMult::Int, oldMult::Int, weights)

Adaption of Algorithm A5.8 from 'The NURBS Book' p. 185.

Modifies solely 'cPts' and 'weights'.
"""
function trimControlPoints!(cPts, kVecOri, degree::Int, pos::Int, pointToRemove::Real, limitedMult::Int, oldMult::Int, weights)

    # --- computation including weights
    if !isempty(weights)
        cPts .*= weights
    else
        weights = ones(size(cPts)) # suboptimal solution: non-existent weights are set to 1
    end

    first = pos - degree
    last  = pos - oldMult

    temp = similar(cPts, 2 * degree + 1)
    wemp = similar(weights, 2 * degree + 1)

    removed = 0

    for t in 0:(limitedMult - 1) # remove limitedMult times the knot

        off = first - 1

        temp[1] = cPts[off]
        temp[last + 2 - off] = cPts[last + 1]

        wemp[1] = weights[off]
        wemp[last + 2 - off] = weights[last + 1]

        i = first
        j = last

        ii = 1
        jj = last - off

        removable = false
        removableW = false

        # --- compute new control points for one removal step
        while j - i > t

            αi = (pointToRemove - kVecOri[i]) / (kVecOri[i + degree + 1 + t] - kVecOri[i])
            αj = (pointToRemove - kVecOri[j - t]) / (kVecOri[j + degree + 1] - kVecOri[j - t])

            temp[ii + 1] = (cPts[i] - (1 - αi) * temp[ii]) / αi
            temp[jj + 1] = (cPts[j] - αj * temp[jj + 2]) / (1 - αj)

            wemp[ii + 1] = (weights[i] - (1 - αi) * wemp[ii]) / αi
            wemp[jj + 1] = (weights[j] - αj * wemp[jj + 2]) / (1 - αj)

            i += 1
            j -= 1

            ii += 1
            jj -= 1
        end

        # --- can knot be removed? Two different cases to be covered
        if j - i < t

            temp[ii] ≈ temp[jj + 2] && (removable = true)
            wemp[ii] ≈ wemp[jj + 2] && (removableW = true)
        else

            αi = (pointToRemove - kVecOri[i]) / (kVecOri[i + degree + 1 + t] - kVecOri[i])
            cPts[i] ≈ (αi * temp[ii + t + 2] + (1 - αi) * temp[ii]) && (removable = true)
            weights[i] ≈ (αi * wemp[ii + t + 2] + (1 - αi) * wemp[ii]) && (removableW = true)
        end

        # --- insert new control points 
        if removable && removableW

            removed += 1

            i = first
            j = last

            while j - i > t

                cPts[i] = temp[i - off + 1]
                cPts[j] = temp[j - off + 1]

                weights[i] = wemp[i - off + 1]
                weights[j] = wemp[j - off + 1]

                i += 1
                j -= 1
            end

        else
            @warn "The knot is only $t times removable."
        end

        first -= 1
        last  += 1
    end

    center  = floor(Int, (2 * pos - degree - oldMult) / 2)
    spStart = center - floor(Int, (removed - 1) / 2)
    spStop  = center + ceil(Int, (removed - 1) / 2)

    # --- remove control points
    splice!(cPts, spStart:spStop)
    splice!(weights, spStart:spStop)

    cPts ./= weights

    return removed # actual number of removed knots
end
