
"""
    numBasisFunctions(basis::Basis)

The number of basis functions is fixed by the knot vector and the degree.

Assumption: the first and last knot vector entry has mulitplicity degree + 1.
"""
numBasisFunctions(basis::Basis) = length(basis.knotVec) - basis.degree - 1


"""
    isValidKnotVector!(kVec)

Check whether the knot vector has only entries in [0, 1] and is in ascending order.

If not the knot vector is modified.
"""
function isValidKnotVector!(kVec)

    # is the vector sorted?
    issorted(kVec) || error("The knot vector is not in ascending order.")

    # is the first element 0?
    kVec[1] == 0.0 || error("The knot vector has to start at 0.")

    # normalize knot vector entries to [0, 1]
    if kVec[end] != 1.0
        kVec ./= maximum(kVec)
        @info "The knot vector is being modified (normalized)."
    end

    return nothing
end


"""
    generateKnotVec(b::Int, degree::Int)

Convenience function to generate a knot vector for 'b' basis functions and a certain 'degree': 

The first and last entry are repeated 'degree'+1 times. Normalized to [0, 1].
"""
function generateKnotVec(b::Int, degree::Int)

    mpc = degree + 1     # required mulitplicity of first and last knot
    m = b + degree + 1   # length of knot vector

    aa = zeros(Int, mpc, 1)
    bb = ones(Int, mpc, 1) * (m - 2 * mpc + 1)

    kVec = Vector{Float64}(undef, m)
    kVec[1:mpc] .= aa
    kVec[(end - mpc + 1):end] .= bb
    kVec[(mpc + 1):(end - mpc)] .= collect(1:(bb[1] - 1))

    kVec ./= maximum(kVec)

    return kVec
end


"""
    Jacobian(Patch::Surface, uEvalpoints, vEvalpoints)

Compute the Jacobian matrix and its (generalized) determinant at the parametric points 'uEvalpoints' and 'vEvalpoints'.

Return
    - J     2-dimensional vector: first for the derivative w.r.t 'u', second w.r.t 'v'
                each vector entry contains a matrix of size (uEvalpoints, vEvalpoints)
                each entry of the matrix is an SVector with the derivatives: SVector(∂x/∂u, ∂y/∂u, ∂y/∂u)

    - dJ    matrix of size (uEvalpoints, vEvalpoints) where each entry is the Jacobi determinant evaluated at the points 'u' and 'v'.

Note: surface points are evaluated but thrown away: maybe change this/make use of it.
"""
function Jacobian(Patch::Surface, uEvalpoints, vEvalpoints)

    S = surfaceDerivativesPoints(Patch, uEvalpoints, vEvalpoints, 1) # first derivatives

    Ju = S[2, 1] # first derivative along u
    Jv = S[1, 2] # first derivative along v

    J = [Ju, Jv] # Jacobi matrix

    # --- Jacobi determinant / surface element
    dJ = Matrix{eltype(Ju[1])}(undef, size(Ju))

    # loop over all points
    for i in axes(Ju, 1)
        for j in axes(Ju, 2)

            xu = Ju[i, j][1]
            yu = Ju[i, j][2]
            zu = Ju[i, j][3]

            xv = Jv[i, j][1]
            yv = Jv[i, j][2]
            zv = Jv[i, j][3]

            dJ[i, j] = sqrt((yu * zv - zu * yv)^2 + (zu * xv - xu * zv)^2 + (xu * yv - yu * xv)^2)
        end
    end

    return J, dJ
end


"""
    spanRanges(Bspl::Bspline, points)

Determine the ranges of the points which lie in each span of the B-spline (assuming normalized open knot vectors).

Return a vector of ranges (one entry per span).
"""
function spanRanges(Bspl::Bspline, points; emptyRanges=false)

    numBasis = numBasisFunctions(Bspl)
    knotSpan = NURBS.findSpan(numBasis, points, Bspl.knotVec) # find for each point the span index
    p = Bspl.degree

    # open knot vector: set first p ranges to 0:-1 (empty range)
    if emptyRanges
        ranges = [0:-1]
        for i in 1:(p - 1)
            push!(ranges, 0:-1)
        end
    else
        ranges = []
    end

    # remaining ranges
    en = 0
    for i in (p + 1):numBasis

        SS = sum(knotSpan .== i) # number of points in span i

        # case: no point lies in the span
        if SS == 0
            push!(ranges, 0:-1) # insert empty range
            continue
        end

        st = en
        en = st + SS

        push!(ranges, (st + 1):en)
    end

    # open knot vector: append p empty ranges
    for i in 1:p
        emptyRanges && push!(ranges, 0:-1)
    end

    return ranges
end
