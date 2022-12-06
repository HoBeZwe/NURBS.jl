
abstract type Basis end

"""
    Bspline{F} <: Basis

B-spline basis.
"""
struct Bspline{F} <: Basis
    degree::Int
    knotVec::Vector{F}
end


"""
    NURB{F} <: Basis

NURBS basis.
"""
struct NURB{F} <: Basis
    degree::Int
    knotVec::Vector{F}
    weights::Vector{F}
end



"""
    Curve{F}

General curve defined by a basis and the control points.
"""
struct Curve{F}
    basis::Basis
    controlPoints::Matrix{SVector{3,F}}
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
