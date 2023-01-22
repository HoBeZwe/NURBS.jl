using NURBS
using Test

using JuliaFormatter
using StaticArrays
using PlotlyJS

# --- testsets
@testset "Testing NURBS functionality" begin

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
        include("knotInsertion.jl")
    end

    @testset "Utils" begin
        include("utils.jl")
    end

    @testset "Test formatting of files" begin
        pkgpath = pkgdir(NURBS)   # path of this package including name
        @test format(pkgpath, overwrite=false)  # check whether files are formatted according to the .JuliaFormatter.toml 
    end
end
