
"""

Remove a knot from a curve 'multiplicity' times.
"""
function removeKnot(C::Curve, pointToRemove::Real, multiplicity::Int)

    knotVecNew, controlPointsNew, weightsNew = removeKnot(
        C.basis.knotVec, C.controlPoints, C.basis.degree, pointToRemove, multiplicity, weights(C.basis)
    )

    return similarCurve(C, C.basis.degree, knotVecNew, controlPointsNew, weightsNew)
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
    actualRemoved = trimControlPoints!(controlPoints, knotVec, degree, oldMultIndices.stop, pointToRemove, limitedMult, oldMult, weights)

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


Adaption of Algorithm A5.8 from 'The NURBS Book' p. 185.

Modifies solely 'cPts' and 'weights'.
"""
function trimControlPoints!(cPts, kVecOri, degree::Int, pos::Int, pointToRemove::Real, limitedMult::Int, oldMult::Int, weights)

    # --- computation including weights
    if !isempty(weights)

        error("Weights not done yet.")

        return nothing
    end

    first = pos - degree
    last  = pos - oldMult

    temp = similar(cPts, 2 * degree + 1)

    removed = 0

    for t in 0:limitedMult-1 # remove limitedMult times the knot

        off = first - 1
 
        temp[1] = cPts[off]
        temp[last + 2 - off] = cPts[last + 1]

        i = first
        j = last

        ii = 1
        jj = last - off

        removable = false

        # --- compute new control points for one removal step
        while j - i > t

            αi = (pointToRemove - kVecOri[i]) / (kVecOri[i + degree + 1 + t] - kVecOri[i])
            αj = (pointToRemove - kVecOri[j - t]) / (kVecOri[j + degree + 1] - kVecOri[j - t])

            temp[ii + 1] = (cPts[i] - (1 - αi) * temp[ii]) / αi
            temp[jj + 1] = (cPts[j] - αj * temp[jj + 2]) / (1 - αj)

            i += 1
            j -= 1

            ii += 1
            jj -= 1
        end

        # --- can knot be removed? Two different cases to be covered
        if j - i < t
            
            temp[ii] ≈ temp[jj + 2] && (removable = true)
        else

            αi = (pointToRemove - kVecOri[i]) / (kVecOri[i + degree + 1 + t] - kVecOri[i])
            cPts[i] ≈ (αi * temp[ii + t + 2] + (1 - αi) * temp[ii]) && (removable = true)
        end

        # --- insert new control points 
        if removable

            removed += 1

            i = first
            j = last

            while j - i > t
                
                cPts[i] = temp[i - off+1]
                cPts[j] = temp[j - off+1]

                i += 1
                j -= 1
            end

        else
            @warn "The knot is only $t times removable."
        end

        first -= 1
        last  += 1
    end    

    center = floor(Int, (2*pos-degree-oldMult)/2 )
    spStart = center - floor(Int, (removed-1) / 2)
    spStop  = center + ceil(Int, (removed-1) / 2)  

    # --- remove control points
    splice!(cPts, spStart:spStop)

    return removed # actual number of removed knots
end
