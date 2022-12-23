
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
