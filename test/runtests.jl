using NURBS
using Test

using JuliaFormatter
using StaticArrays
using LinearAlgebra
using PlotlyJS
using FileIO

# --- testsets
@testset verbose = true "Testing NURBS functionality" begin

    @testset "Bases" begin
        include("bases.jl")
    end

    @testset "Curves" begin
        include("curves_Bspline.jl")
        include("curves_nurbs.jl")
    end

    @testset "Surfaces" begin
        include("surfaces_Bspline.jl")
        include("surfaces_nurbs.jl")
    end

    @testset "Fundamental Operations" begin
        include("fundamentalOperations/knotInsertion.jl")
        include("fundamentalOperations/splitting.jl")
        include("fundamentalOperations/knotRemoval.jl")
    end

    @testset "Utils" begin
        include("utils.jl")
    end

    @testset "Connectivity" begin
        include("connectivity/interfaces.jl")
        include("connectivity/bezierMesh.jl")
    end

    @testset "Test formatting of files" begin
        pkgpath = pkgdir(NURBS)   # path of this package including name
        @test format(pkgpath, overwrite=false)  # check whether files are formatted according to the .JuliaFormatter.toml 
    end
end
