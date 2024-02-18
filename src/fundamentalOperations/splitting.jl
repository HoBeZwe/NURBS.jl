
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

    @suppress begin # suppress info output when constructing the curves

        for (i, splitPoint) in enumerate(splits)

            checkRange(splitPoint) # is input a valid one?

            kVec1, kVec, w1, w, cPts1, cPts = splitData(p, w, kVec, cPts, splitPoint)

            push!(cVec, similarCurve(C, p, kVec1, cPts1, w1))
        end

        push!(cVec, similarCurve(C, p, kVec, cPts, w))

    end

    return cVec
end


"""
    Base.split(C::Curve, splitPoint=0.5)

Split a curve at a single parametric point in the range ]0,1[ by inserting a single knot multiple times.
"""
function Base.split(C::Curve, splitPoint=0.5)

    checkRange(splitPoint) # is input a valid one?

    p = C.basis.degree
    w = weights(C.basis)
    kVec = C.basis.knotVec
    cPts = C.controlPoints

    kVec1, kVec2, w1, w2, cPts1, cPts2 = splitData(p, w, kVec, cPts, splitPoint)

    C1 = similarCurve(C, C.basis.degree, kVec1, cPts1, w1)
    C2 = similarCurve(C, C.basis.degree, kVec2, cPts2, w2)

    return C1, C2
end


"""
    Base.split(C::Curve, n::Int)

Split a curve in the parametric domain equally to obtain n curves.
"""
function Base.split(S::Surface, n::Int, m::Int)

    sep = 1.0 / n
    splitsU = collect(sep:sep:(1.0 - sep))

    sep = 1.0 / m
    splitsV = collect(sep:sep:(1.0 - sep))

    return split(S; U=splitsU, V=splitsV)
end


"""
    Base.split(S::SurfaceT; U=[], V=[]) where {SurfaceT<:Surface}

Split a surface at multiple parametric points in U and V, each in the range ]0,1[, by inserting knots.
"""
function Base.split(S::SurfaceT; U=[], V=[]) where {SurfaceT<:Surface}

    issorted(U) || error("The parametric points have to be sorted.")
    issorted(V) || error("The parametric points have to be sorted.")
    allunique(U) || error("There are split points occuring more than once.")
    allunique(V) || error("There are split points occuring more than once.")

    sVec = SurfaceT[]

    @suppress begin # suppress info output when constructing the surfaces

        # --- all splits along u
        Su = splitU(S, U)

        # --- split each surface along v
        for (i, su) in enumerate(Su)
            append!(sVec, splitV(su, V))
        end
    end

    return sVec
end


"""
    splitU(S::SurfaceT, splits::Vector) where {SurfaceT<:Surface}

Split a surface along the u-direction at multiple parametric points, each in the range ]0,1[, by inserting knots.

The splits vector needs to be sorted and each value should only occur once (both are not checked since this function is meant for internal use).
"""
function splitU(S::SurfaceT, splits::Vector) where {SurfaceT<:Surface}

    sVec = SurfaceT[]

    degree = S.uBasis.degree
    kVec = S.uBasis.knotVec
    cPts = S.controlPoints
    wAux = weights(S)

    for (i, splitPoint) in enumerate(splits)

        checkRange(splitPoint) # is input a valid one?

        # --- split knot vector: increase multiplicity of split entry to degree
        kVec1, kVec2, oldMultIndices, oldMult, limitedMult = splitKnots(degree, kVec, splitPoint)

        # --- split control points + weights
        N1 = numBasisFunctions(kVec1, degree)
        N2 = numBasisFunctions(kVec2, degree)
        M = size(cPts, 2)

        wMat1 = similar(weights(S), N1, M)
        wMat2 = similar(weights(S), N2, M)

        cPts1 = similar(cPts, N1, M)
        cPts2 = similar(cPts, N2, M)

        for i in eachindex(cPts[1, :])

            cPts1[:, i], cPts2[:, i], wMod = splitControlPoints(
                degree, weights(wAux, :, i), cPts[:, i], kVec, oldMultIndices, oldMult, limitedMult, splitPoint
            )

            if !isempty(wAux)
                wMat1[:, i], wMat2[:, i] = splitWeights(wMod, oldMultIndices.start)
            end
        end

        cPts = cPts2
        kVec = kVec2
        isempty(wAux) || (wAux = wMat2)

        push!(sVec, similarSurface(S, degree, S.vBasis.degree, kVec1, S.vBasis.knotVec, cPts1, wMat1))
    end

    push!(sVec, similarSurface(S, degree, S.vBasis.degree, kVec, S.vBasis.knotVec, cPts, wAux))

    return sVec
end


"""
    splitV(S::SurfaceT, splits::Vector) where {SurfaceT<:Surface}

Split a surface along the v-direction at multiple parametric points, each in the range ]0,1[, by inserting knots.

The splits vector needs to be sorted and each value should only occur once (both are not checked since this function is meant for internal use).
"""
function splitV(S::SurfaceT, splits::Vector) where {SurfaceT<:Surface}

    sVec = SurfaceT[]

    degree = S.vBasis.degree
    kVec = S.vBasis.knotVec
    cPts = S.controlPoints
    wAux = weights(S)

    for (i, splitPoint) in enumerate(splits)

        checkRange(splitPoint) # is input a valid one?

        # --- split knot vector: increase multiplicity of split entry to degree
        kVec1, kVec2, oldMultIndices, oldMult, limitedMult = splitKnots(degree, kVec, splitPoint)

        # --- split control points + weights
        N1 = numBasisFunctions(kVec1, degree)
        N2 = numBasisFunctions(kVec2, degree)
        M = size(cPts, 1)

        wMat1 = similar(weights(S), M, N1)
        wMat2 = similar(weights(S), M, N2)

        cPts1 = similar(cPts, M, N1)
        cPts2 = similar(cPts, M, N2)

        for i in eachindex(cPts[:, 1])

            cPts1[i, :], cPts2[i, :], wMod = splitControlPoints(
                degree, weights(wAux, i, :), cPts[i, :], kVec, oldMultIndices, oldMult, limitedMult, splitPoint
            )

            if !isempty(wAux)
                wMat1[i, :], wMat2[i, :] = splitWeights(wMod, oldMultIndices.start)
            end
        end

        cPts = cPts2
        kVec = kVec2
        isempty(wAux) || (wAux = wMat2)

        push!(sVec, similarSurface(S, S.uBasis.degree, degree, S.uBasis.knotVec, kVec1, cPts1, wMat1))
    end

    push!(sVec, similarSurface(S, S.uBasis.degree, degree, S.uBasis.knotVec, kVec, cPts, wAux))

    return sVec
end


"""
    splitData(degree::Int, wOri, kVecOri, ctrlPtsOri, splitPoint=0.5) 

Split the underlying data of a curve at a single parametric point in the range ]0,1[ by inserting a single knot multiple times.
"""
function splitData(degree::Int, wOri, kVecOri, ctrlPtsOri, splitPoint)

    # --- split knot vector: increase multiplicity of split entry to degree
    kVec1, kVec2, oldMultIndices, oldMult, limitedMult = splitKnots(degree, kVecOri, splitPoint)

    # --- split control points
    cPts1, cPts2, wMod = splitControlPoints(degree, wOri, ctrlPtsOri, kVecOri, oldMultIndices, oldMult, limitedMult, splitPoint)

    # --- split weights
    w1, w2 = splitWeights(wMod, oldMultIndices.start)

    return kVec1, kVec2, w1, w2, cPts1, cPts2
end


"""
    splitControlPoints(degree::Int, wOri, ctrlPtsOri, kVecOri, oldMultIndices, oldMult, limitedMult, splitPoint)

Split the control points (and the weights) by inserting a knot multiple times.
"""
function splitControlPoints(degree::Int, wOri, ctrlPtsOri, kVecOri, oldMultIndices, oldMult, limitedMult, splitPoint)

    # --- copy to modify
    wMod = deepcopy(wOri)
    ctrlPtsMod = deepcopy(ctrlPtsOri)

    # --- insert knot
    extendControlPoints!(ctrlPtsMod, kVecOri, degree, oldMultIndices.stop, splitPoint, limitedMult, oldMult, wMod)

    # --- split
    cPts1 = ctrlPtsMod[1:(oldMultIndices.start - 1)]
    cPts2 = ctrlPtsMod[(oldMultIndices.start - 1):end]

    return cPts1, cPts2, wMod
end


"""
    splitKnots(degree::Int, kVecOri, splitPoint)

Split knot vector according to the split point.
"""
function splitKnots(degree::Int, kVecOri, splitPoint)

    # --- copy to modify
    knotVecNew = deepcopy(kVecOri)

    # --- insert knot
    oldMultIndices, oldMult, limitedMult = extendKnotVector!(knotVecNew, degree, splitPoint, degree)

    # --- split
    kVec1 = [knotVecNew[1:(oldMultIndices.start + degree - 1)]; splitPoint]
    kVec2 = [splitPoint; knotVecNew[(oldMultIndices.start):end]]

    return kVec1, kVec2, oldMultIndices, oldMult, limitedMult
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
