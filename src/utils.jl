
"""
    numBasisFunctions(basis::Basis)

The number of basis functions is fixed by the knot vector and the degree.

Assumption: the first and last knot vector entry has mulitplicity degree + 1.
"""
numBasisFunctions(basis::Basis) = length(basis.knotVec) - basis.degree - 1


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
    jacobian(Patch::Surface, uEvalpoints, vEvalpoints)

Compute the Jacobian matrix and its (generalized) determinant at the parametric points 'uEvalpoints' and 'vEvalpoints'.

Return
    - J     2-dimensional vector: first for the derivative w.r.t 'u', second w.r.t 'v'
                each vector entry contains a matrix of size (uEvalpoints, vEvalpoints)
                each entry of the matrix is an SVector with the derivatives: SVector(∂x/∂u, ∂y/∂u, ∂y/∂u)

    - dJ    matrix of size (uEvalpoints, vEvalpoints) where each entry is the Jacobi determinant evaluated at the points 'u' and 'v'.

Note: surface points are evaluated but thrown away: maybe change this/make use of it.
"""
function jacobian(Patch::Surface, uEvalpoints, vEvalpoints)

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
