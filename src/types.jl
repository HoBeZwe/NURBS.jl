
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
    CurrySchoenberg{F} <: Basis{F}

Normalized B-spline basis.
"""
struct CurrySchoenberg{F} <: Basis{F}
    degree::Int
    knotVec::Vector{F}

    # inner constructor: check whether provided knot vector a valid one
    function CurrySchoenberg(degree::Int, knotVec::Vector{F}) where {F}
        isValidKnotVector!(knotVec)
        new{F}(degree, knotVec)
    end
end


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



abstract type Shape{F} end

abstract type Curve{F} <: Shape{F} end

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




abstract type Surface{F} <: Shape{F} end

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





"""
    similarCurve(curve::BsplineCurve, p::Int, kVec, cPts, w)

Construct B-spline curve from underlying data: ignore empty weights.
"""
function similarCurve(curve::BsplineCurve, p::Int, kVec, cPts, w)

    return BsplineCurve(Bspline(p, kVec), cPts)
end


"""
    similarCurve(curve::NURBScurve, p::Int, kVec, cPts, w)

Construct NURBS curve from underlying data.
"""
function similarCurve(curve::NURBScurve, p::Int, kVec, cPts, w)

    return NURBScurve(NURB(p, kVec, w), cPts)
end


"""
    similarSurface(surface::BsplineSurface, p::Int, uVec, vVec, cPts, w)

Construct B-spline surface from underlying data: ignore empty weights.
"""
function similarSurface(surface::BsplineSurface, p::Int, uVec, vVec, cPts, w)

    return BsplineSurface(Bspline(p, uVec), Bspline(p, vVec), cPts)
end


"""
    similarSurface(surface::NURBSsurface, p::Int, uVec, vVec, cPts, weights)

Construct NURBS surface from underlying data.
"""
function similarSurface(surface::NURBSsurface, p::Int, uVec, vVec, cPts, weights)

    return NURBSsurface(Bspline(p, uVec), Bspline(p, vVec), cPts, weights)
end
