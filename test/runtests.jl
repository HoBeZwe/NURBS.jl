
using TestItemRunner

@testitem "Bases" begin
    include("bases.jl")
end

@testitem "Curves" begin
    include("curves_Bspline.jl")
    include("curves_nurbs.jl")
end

@testitem "Surfaces" begin
    include("surfaces_Bspline.jl")
    include("surfaces_nurbs.jl")
end

@testitem "Fundamental Operations" begin
    include("fundamentalOperations/knotInsertion.jl")
    include("fundamentalOperations/splitting.jl")
    include("fundamentalOperations/knotRemoval.jl")
end

@testitem "Utils" begin
    include("utils.jl")
end

@testitem "Connectivity" begin
    include("connectivity/interfaces.jl")
    include("connectivity/bezierMesh.jl")
end

@testitem "Formatting of files" begin
    using JuliaFormatter
    pkgpath = pkgdir(NURBS)                 # path of this package including name
    @test format(pkgpath, overwrite=false)  # check whether files are formatted according to the .JuliaFormatter.toml 
end



@run_package_tests verbose = true
