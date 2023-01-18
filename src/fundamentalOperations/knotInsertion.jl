
"""
    insertKnot(knotVecOrig, controlPointsOrig, p::Int, newParametricPoint::Real, multiplicity::Int)

Insert a new parametric point with a given mulitplicity (insert the point multiple times) for the polynomial degree 'degree'. Return the resutling knot vector and the new control points.

The provided knot vector and control-point vector are NOT modified.
"""
function insertKnot(knotVecOrig, controlPointsOrig, degree::Int, newParametricPoint::Real, multiplicity::Int)

    # --- make copies
    knotVecNew       = deepcopy(knotVecOrig)
    controlPointsNew = deepcopy(controlPointsOrig)

    # --- modify copies
    insertKnot!(knotVecNew, controlPointsNew, degree, newParametricPoint, multiplicity)

    return knotVecNew, controlPointsNew
end


"""
    insertKnot!(knotVecOrig, controlPointsOrig, p::Int, newParametricPoint::Real, multiplicity::Int)

Insert a new parametric point with a given mulitplicity (insert the point multiple times) for the polynomial degree 'degree'. Return the resutling knot vector and the new control points.

The provided knot vector and control-point vector are modified.
"""
function insertKnot!(knotVecOrig, controlPointsOrig, degree::Int, newParametricPoint::Real, multiplicity::Int)

    # --- is provide knot vector a valid one?
    isValidKnotVector!(knotVecOrig)


    # --- modify knot vector
    oldMultIndices, oldMult = extendKnotVector!(knotVecOrig, degree, newParametricPoint, multiplicity)
    
    
    # --- modify control points


    return nothing
end


"""
    extendKnotVector!(knotVecOrig, newParametricPoint::Real, multiplicity::Int)

Insert the new parametric point with a given mulitplicity (insert the point multiple times).
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
