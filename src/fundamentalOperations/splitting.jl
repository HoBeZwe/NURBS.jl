
"""
    Base.split(C::Curve, n::Int)

Split a curve in the parametric domain equally to obtain n curves.
"""
function Base.split(C::Curve, n::Int)

    sep = 1.0 / n
    splits = collect(sep:sep:(1.0 - sep))

    return split(C, splits)
end



"""
    Base.split(C::CurveT, splits::Vector) where {CurveT<:Curve}

Split a curve at multiple parametric points, each in the range ]0,1[, by inserting knots.

The splits vector needs to be sorted and each value should only occur once.
"""
function Base.split(C::CurveT, splits::Vector) where {CurveT<:Curve}

    issorted(splits) || error("The parametric points have to be sorted.")
    allunique(splits) || error("There are split points occuring more than once.")

    cVec = CurveT[]

    p = C.basis.degree
    w = weights(C.basis)
    kVec = C.basis.knotVec
    cPts = C.controlPoints

    for (i, splitPoint) in enumerate(splits)

        (splitPoint > 0.0 && splitPoint < 1.0) || error("A parametric split point is outside ]0,1[.")

        kVec1, kVec, w1, w, cPts1, cPts = splitData(p, w, kVec, cPts, splitPoint)

        push!(cVec, splittedCurve(C, p, kVec1, cPts1, w1))
    end

    push!(cVec, splittedCurve(C, p, kVec, cPts, w))

    return cVec
end


"""
    Base.split(C::Curve, splitPoint=0.5)

Split a curve at a single parametric point in the range ]0,1[ by inserting a single knot.
"""
function Base.split(C::Curve, splitPoint=0.5)

    (splitPoint > 0.0 && splitPoint < 1.0) || error("The parametric split point is outside ]0,1[.")

    p = C.basis.degree
    w = weights(C.basis)
    kVec = C.basis.knotVec
    cPts = C.controlPoints

    kVec1, kVec2, w1, w2, cPts1, cPts2 = splitData(p, w, kVec, cPts, splitPoint)

    C1 = splittedCurve(C, C.basis.degree, kVec1, cPts1, w1)
    C2 = splittedCurve(C, C.basis.degree, kVec2, cPts2, w2)

    return C1, C2
end


"""
    splitData(pOri::Int, wOri, kVecOri, ctrlPtsOri, splitPoint=0.5) 

Split the underlying data of a curve at a single parametric point in the range ]0,1[ by inserting a single knot.
"""
function splitData(pOri::Int, wOri, kVecOri, ctrlPtsOri, splitPoint)

    # --- copies to modify
    wNew = deepcopy(wOri)
    knotVecNew = deepcopy(kVecOri)
    ctrlPtsNew = deepcopy(ctrlPtsOri)


    # --- split knot vector: increase multiplicity of split entry to pOri
    oldMultIndices, oldMult, limitedMult = extendKnotVector!(knotVecNew, pOri, splitPoint, pOri)

    kVec1 = [knotVecNew[1:(oldMultIndices.start + pOri - 1)]; splitPoint]
    kVec2 = [splitPoint; knotVecNew[(oldMultIndices.start):end]]


    # --- split control points
    extendControlPoints!(ctrlPtsNew, kVecOri, pOri, oldMultIndices.stop, splitPoint, limitedMult, oldMult, wNew)

    cPts1 = ctrlPtsNew[1:(oldMultIndices.start - 1)]
    cPts2 = ctrlPtsNew[(oldMultIndices.start - 1):end]


    # --- split weights
    w1, w2 = splitWeights(wNew, oldMultIndices.start)


    return kVec1, kVec2, w1, w2, cPts1, cPts2
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


"""
    splittedCurve(curve::BsplineCurve, p::Int, kVec, cPts, w)

Construct B-spline curve from underlying data: ignore empty weights.
"""
function splittedCurve(curve::BsplineCurve, p::Int, kVec, cPts, w)

    return BsplineCurve(Bspline(p, kVec), cPts)
end


"""
    splittedCurve(curve::NURBScurve, p::Int, kVec, cPts, w)

Construct NURBS curve from underlying data.
"""
function splittedCurve(curve::NURBScurve, p::Int, kVec, cPts, w)

    return NURBScurve(NURB(p, kVec, w), cPts)
end
