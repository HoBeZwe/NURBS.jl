
"""

Split a curve by inserting knots.
"""
function Base.split(C::CurveT, splits::Vector=[0.5]) where {CurveT<:Curve}

    pOri = C.basis.degree
    wNew = deepcopy(weights(C.basis))
    knotVecNew = deepcopy(C.basis.knotVec)
    ctrlPtsNew = deepcopy(C.controlPoints)


    # --- split knot vector: increase multiplicity of split entry to pOri
    oldMultIndices, oldMult, limitedMult = extendKnotVector!(knotVecNew, pOri, splits[1], pOri)

    kVec1 = [knotVecNew[1:(oldMultIndices.start + pOri - 1)]; splits[1]]
    kVec2 = [splits[1]; knotVecNew[(oldMultIndices.start):end]]


    # --- split control points
    extendControlPoints!(ctrlPtsNew, C.basis.knotVec, pOri, oldMultIndices.stop, splits[1], limitedMult, oldMult, wNew)

    cPts1 = ctrlPtsNew[1:(oldMultIndices.start - 1)]
    cPts2 = ctrlPtsNew[(oldMultIndices.start - 1):end]


    # --- split weights
    w1, w2 = splitWeights(wNew, oldMultIndices.start)


    # --- return curves
    C1 = splittedCurve(C, pOri, kVec1, cPts1, w1)
    C2 = splittedCurve(C, pOri, kVec2, cPts2, w2)

    return C1, C2
end


"""
    splitWeights(weights::Vector{T}, Ind::Int) where {T}

Split the weights including the case when there are no weights.
"""
function splitWeights(weights::Vector{T}, Ind::Int) where {T}

    isempty(weights) && return weights, weights

    w1 = weights[1:(Ind - 1)]
    w2 = weights[(Ind - 1):end]

    return w1, w2
end


function splittedCurve(curve::BsplineCurve, p::Int, kVec, cPts, w)

    return BsplineCurve(Bspline(p,kVec), cPts)
end

function splittedCurve(curve::NURBScurve, p::Int, kVec, cPts, w)

    return NURBScurve(NURB(p, kVec, w), cPts)
end