
abstract type Basis{F} end

"""
    Bspline{F} <: Basis{F}

B-spline basis.
"""
struct Bspline{F} <: Basis{F}
    degree::Int
    knotVec::Vector{F}

    # inner constructor: check whether provided knot vector a valid one
    function Bspline(degree::Int, knotVec::Vector{F}) where {F}
        isValidKnotVector!(knotVec)
        new{F}(degree, knotVec)
    end
end


"""
    BsplineDerivatives{F} <: Basis{F}

B-spline basis including all derivatives up to divMax.
"""
struct BsplineDerivatives{F} <: Basis{F}
    degree::Int
    knotVec::Vector{F}
    divMax::Int
end

derivatives(B::Bspline, k::Int) = BsplineDerivatives(B.degree, B.knotVec, k)
∂(B::Bspline) = BsplineDerivatives(B.degree, B.knotVec, 1) # first derivative


"""
    NURB{F} <: Basis{F}

NURBS basis.
"""
struct NURB{F} <: Basis{F}
    degree::Int
    knotVec::Vector{F}
    weights::Vector{F}

    # inner constructor: perform some checks
    function NURB(degree::Int, knotVec::Vector{F}, weights::Vector{F}) where {F}
        isValidKnotVector!(knotVec)      # is provided knot vector a valid one?

        B = length(knotVec) - degree - 1 # number of basis functions
        length(weights) == B || error("The length of the weights vector has to match the number of basis functions ($B).")

        new{F}(degree, knotVec, weights)
    end
end




abstract type Curve{F} end

"""
    BsplineCurve{F} <: Curve{F}

B-spline curve defined by the basis and the control points.
"""
struct BsplineCurve{F} <: Curve{F}
    basis::Bspline{F}
    controlPoints::Vector{SVector{3,F}}
end

"""
    NURBScurve{F} <: Curve{F}

B-spline curve defined by the basis and the control points.
"""
struct NURBScurve{F} <: Curve{F}
    basis::NURB{F}
    controlPoints::Vector{SVector{3,F}}
end




abstract type Surface{F} end

"""
    BsplineSurface{F} <: Surface{F}

Surface defined by a B-spline basis and the control points.
"""
struct BsplineSurface{F} <: Surface{F}
    uBasis::Bspline{F}
    vBasis::Bspline{F}
    controlPoints::Matrix{SVector{3,F}}     # control points in (u, v)- direction
end


"""
    NURBSsurface{F} <: Surface{F}

Surface defined by a B-spline basis, the control points, and the weights.
"""
struct NURBSsurface{F} <: Surface{F}
    uBasis::Bspline{F}
    vBasis::Bspline{F}
    controlPoints::Matrix{SVector{3,F}}     # control points in (u, v)- direction
    weights::Matrix{F}                      # weights for each control point in same ordering
end
