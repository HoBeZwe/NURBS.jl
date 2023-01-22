
abstract type Basis end

"""
    Bspline{F} <: Basis

B-spline basis.
"""
struct Bspline{F} <: Basis
    degree::Int
    knotVec::Vector{F}

    # inner constructor: check whether provided knot vector a valid one
    function Bspline(degree::Int, knotVec::Vector{F}) where {F}
        isValidKnotVector!(knotVec)
        new{F}(degree, knotVec)
    end
end


"""
    NURB{F} <: Basis

NURBS basis.
"""
struct NURB{F} <: Basis
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




abstract type Curve end

"""
    BsplineCurve{F} <: Curve

B-spline curve defined by the basis and the control points.
"""
struct BsplineCurve{F} <: Curve
    basis::Bspline{F}
    controlPoints::Vector{SVector{3,F}}
end

"""
    NURBScurve{F} <: Curve

B-spline curve defined by the basis and the control points.
"""
struct NURBScurve{F} <: Curve
    basis::NURB{F}
    controlPoints::Vector{SVector{3,F}}
end




abstract type Surface end

"""
    BsplineSurface{F} <: Surface

Surface defined by a B-spline basis and the control points.
"""
struct BsplineSurface{F} <: Surface
    uBasis::Bspline{F}
    vBasis::Bspline{F}
    controlPoints::Matrix{SVector{3,F}}     # control points in (u, v)- direction
end


"""
    NURBSsurface{F} <: Surface

Surface defined by a B-spline basis, the control points, and the weights.
"""
struct NURBSsurface{F} <: Surface
    uBasis::Bspline{F}
    vBasis::Bspline{F}
    controlPoints::Matrix{SVector{3,F}}     # control points in (u, v)- direction
    weights::Matrix{F}                      # weights for each control point in same ordering
end
