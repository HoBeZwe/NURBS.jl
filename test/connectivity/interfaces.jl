
using StaticArrays
using LinearAlgebra
using PlotlyJS
using FileIO

@testset "Interfaces" begin

    @testset "Sphere" begin

        Patches = readMultipatch("assets/sphere.dat")

        interfaces, commonVtxs = identifyInterfaces(Patches)

        # manually verified values for the sphere
        p1vec = [1 1 1 1 2 2 2 3 3 3 4 4]
        p2vec = [2 4 5 6 3 5 6 4 5 6 5 6]
        l1vec = [4 2 3 1 4 3 1 3 2 4 2 4]
        l2vec = [2 3 2 3 1 1 4 1 4 1 3 2]
        revec = [false true true false true false true false false true true false]
        orvec = [1 1 1 1 1 1 1 1 1 1 1 1]

        for (ind, ifc) in enumerate(interfaces)

            p1, p2 = patchIDs(ifc)
            l1, l2 = localEdges(ifc)

            @test p1 == p1vec[ind]
            @test p2 == p2vec[ind]
            @test l1 == l1vec[ind]
            @test l2 == l2vec[ind]

            @test ifc.reverse == revec[ind]
            @test ifc.orientation == orvec[ind]
        end

        @test isempty(commonVtxs) # there are no patches that share only one vertex
    end


    @testset "Patches with corner touch" begin

        Patches = readMultipatch("assets/fourPatches.dat")

        interfaces, commonVtxs = identifyInterfaces(Patches)

        # manually verified values for the sphere
        p1vec = [1 2 2]
        p2vec = [4 3 4]
        l1vec = [3 2 2]
        l2vec = [3 4 3]

        for (ind, cvtx) in enumerate(commonVtxs)

            @test patchID(cvtx[1]) == p1vec[ind]
            @test patchID(cvtx[2]) == p2vec[ind]

            @test localEdge(cvtx[1]) == l1vec[ind]
            @test localEdge(cvtx[2]) == l2vec[ind]
        end
    end


    @testset "Not supported case" begin

        Patches = readMultipatch("assets/sphere.dat")
        Patches = [Patches[1], Patches[1]] # two identical patches

        @test_throws ErrorException(
            "Patches touch at more than two corner points (4 corner points). So far only one common edge between patches is supported!\n",
        ) identifyInterfaces(Patches)

        @test_throws ErrorException("Point not found!") NURBS.getGlobalIndex(SVector(2.1, 0, 0), [SVector(2.4, 1.1, 0.3)])
    end
end
