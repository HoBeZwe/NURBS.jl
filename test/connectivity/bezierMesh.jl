
using StaticArrays
using LinearAlgebra
using PlotlyJS
using FileIO

@testset "Bezier mesh" begin

    @testset "Sphere" begin

        patch1Le = [
            [Set([48]), Set([36]), Set([2]), Set([4])],
            [Set([1]), Set([33]), Set([3]), Set([5])],
            [Set([2]), Set([30]), Set([39]), Set([6])],
            [Set([51]), Set([1]), Set([5]), Set([7])],
            [Set([4]), Set([2]), Set([6]), Set([8])],
            [Set([5]), Set([3]), Set([38]), Set([9])],
            [Set([54]), Set([4]), Set([8]), Set([10])],
            [Set([7]), Set([5]), Set([9]), Set([11])],
            [Set([8]), Set([6]), Set([37]), Set([12])],
        ]

        patch2Le = [
            [Set([54]), Set([7]), Set([11]), Set([13])],
            [Set([10]), Set([8]), Set([12]), Set([14])],
            [Set([11]), Set([9]), Set([37]), Set([15])],
            [Set([53]), Set([10]), Set([14]), Set([16])],
            [Set([13]), Set([11]), Set([15]), Set([17])],
            [Set([14]), Set([12]), Set([40]), Set([18])],
            [Set([52]), Set([13]), Set([17]), Set([25])],
            [Set([16]), Set([14]), Set([18]), Set([22])],
            [Set([17]), Set([15]), Set([43]), Set([19])],
        ]

        patch5Le = [
            [Set([12]), Set([9]), Set([38]), Set([40])],
            [Set([37]), Set([6]), Set([39]), Set([41])],
            [Set([38]), Set([3]), Set([30]), Set([42])],
            [Set([15]), Set([37]), Set([41]), Set([43])],
            [Set([40]), Set([38]), Set([42]), Set([44])],
            [Set([41]), Set([39]), Set([29]), Set([45])],
            [Set([18]), Set([40]), Set([44]), Set([19])],
            [Set([43]), Set([41]), Set([45]), Set([20])],
            [Set([44]), Set([42]), Set([28]), Set([21])],
        ]



        patch1Lv = [
            [Set([]), Set([33]), Set([5]), Set([51])],
            [Set([36]), Set([30]), Set([6]), Set([4])],
            [Set([33]), Set([]), Set([38]), Set([5])],
            [Set([48]), Set([2]), Set([8]), Set([54])],
            [Set([1]), Set([3]), Set([9]), Set([7])],
            [Set([2]), Set([39]), Set([37]), Set([8])],
            [Set([51]), Set([5]), Set([11]), Set([])],
            [Set([4]), Set([6]), Set([12]), Set([10])],
            [Set([5]), Set([38]), Set([]), Set([11])],
        ]

        patch2Lv = [
            [Set([]), Set([8]), Set([14]), Set([53])],
            [Set([7]), Set([9]), Set([15]), Set([13])],
            [Set([8]), Set([]), Set([40]), Set([14])],
            [Set([54]), Set([11]), Set([17]), Set([52])],
            [Set([10]), Set([12]), Set([18]), Set([16])],
            [Set([11]), Set([37]), Set([43]), Set([17])],
            [Set([53]), Set([14]), Set([22]), Set([])],
            [Set([13]), Set([15]), Set([19]), Set([25])],
            [Set([14]), Set([40]), Set([]), Set([22])],
        ]

        patch3Lv = [
            [Set([]), Set([44]), Set([23]), Set([17])],
            [Set([43]), Set([45]), Set([24]), Set([22])],
            [Set([44]), Set([]), Set([31]), Set([23])],
            [Set([18]), Set([20]), Set([26]), Set([16])],
            [Set([19]), Set([21]), Set([27]), Set([25])],
            [Set([20]), Set([28]), Set([34]), Set([26])],
            [Set([17]), Set([23]), Set([49]), Set([])],
            [Set([22]), Set([24]), Set([46]), Set([52])],
            [Set([23]), Set([31]), Set([]), Set([49])],
        ]

        patch4Lv = [
            [Set([]), Set([42]), Set([32]), Set([24])],
            [Set([45]), Set([39]), Set([33]), Set([31])],
            [Set([42]), Set([]), Set([2]), Set([32])],
            [Set([21]), Set([29]), Set([35]), Set([27])],
            [Set([28]), Set([30]), Set([36]), Set([34])],
            [Set([29]), Set([3]), Set([1]), Set([35])],
            [Set([24]), Set([32]), Set([47]), Set([])],
            [Set([31]), Set([33]), Set([48]), Set([46])],
            [Set([32]), Set([2]), Set([]), Set([47])],
        ]

        patch5Lv = [
            [Set([]), Set([6]), Set([41]), Set([15])],
            [Set([9]), Set([3]), Set([42]), Set([40])],
            [Set([6]), Set([]), Set([29]), Set([41])],
            [Set([12]), Set([38]), Set([44]), Set([18])],
            [Set([37]), Set([39]), Set([45]), Set([43])],
            [Set([38]), Set([30]), Set([28]), Set([44])],
            [Set([15]), Set([41]), Set([20]), Set([])],
            [Set([40]), Set([42]), Set([21]), Set([19])],
            [Set([41]), Set([29]), Set([]), Set([20])],
        ]

        patch6Lv = [
            [Set([]), Set([35]), Set([50]), Set([26])],
            [Set([34]), Set([36]), Set([51]), Set([49])],
            [Set([35]), Set([]), Set([4]), Set([50])],
            [Set([27]), Set([47]), Set([53]), Set([25])],
            [Set([46]), Set([48]), Set([54]), Set([52])],
            [Set([47]), Set([1]), Set([7]), Set([53])],
            [Set([26]), Set([50]), Set([13]), Set([])],
            [Set([49]), Set([51]), Set([10]), Set([16])],
            [Set([50]), Set([4]), Set([]), Set([13])],
        ]

        localVertsRef = [patch1Lv; patch2Lv; patch3Lv; patch4Lv; patch5Lv; patch6Lv]


        Patches = readMultipatch("assets/sphere.dat")

        cty = Connectivity(Patches, 3, 3) # construct connectivity
        badj = cty.bezierAdjacency

        # vertex adjacency
        for i in eachindex(badj) # compare with manually created table
            @test badj[i].atLocalVerts == localVertsRef[i]
        end

        # edge adjacency
        for i in 1:9 # compare with manually created table for patch 1
            @test badj[i].atLocalEdges == patch1Le[i]
        end
        for i in 10:18 # compare with manually created table for patch 2
            @test badj[i].atLocalEdges == patch2Le[i - 9]
        end
        for i in 37:45 # compare with manually created table for patch 5
            @test badj[i].atLocalEdges == patch5Le[i - 36]
        end
    end


    @testset "4 Patches with corner touch" begin

        LocalV1to8 = [
            [Set([]), Set([]), Set([6]), Set([])],
            [Set([]), Set([]), Set([7]), Set([5])],
            [Set([]), Set([]), Set([8]), Set([6])],
            [Set([]), Set([]), Set([37]), Set([7])],
            [Set([]), Set([2]), Set([10]), Set([])],
            [Set([1]), Set([3]), Set([11]), Set([9])],
            [Set([2]), Set([4]), Set([12]), Set([10])],
            [Set([3]), Set([33]), Set([41]), Set([11])],
        ]

        LocalV13to16 = [
            [Set([]), Set([10]), Set([18]), Set([])],
            [Set([9]), Set([11]), Set([19]), Set([17])],
            [Set([10]), Set([12]), Set([20]), Set([18])],
            [Set([11]), Set([41]), Set([64]), Set([19])],
        ]

        LocalV20 = [Set([15]), Set([64, 45]), Set([]), Set([23])]
        LocalV45 = [Set([12]), Set([42]), Set([63]), Set([20])]
        LocalV64 = [Set([59]), Set([]), Set([20, 16]), Set([46])]



        Patches = readMultipatch("assets/fourPatches.dat")

        cty = Connectivity(Patches, 4, 4) # construct connectivity
        badj = cty.bezierAdjacency

        # vertex adjacency
        for i in 1:8 # compare with manually created table
            @test badj[i].atLocalVerts == LocalV1to8[i]
        end
        for i in 13:16 # compare with manually created table
            @test badj[i].atLocalVerts == LocalV13to16[i - 12]
        end
        @test badj[20].atLocalVerts == LocalV20
        @test badj[45].atLocalVerts == LocalV45
        @test badj[64].atLocalVerts == LocalV64
    end

    @testset "Corner cases" begin

        @test NURBS.touchingEdge(1, 1, 17, false, 4, 4) == -1
        @test NURBS.touchingVert(1, 1, 17, false, 4, 4) == (-1, -1)
        @test NURBS.cornerBezierCell(17, 1, 1, 1) == -1

        pInd, cellU, cellV = NURBS.cellLin2Cart(17, 4, 4)
        @test NURBS.cellCart2Lin(cellU, cellV, pInd, 4, 4) == 17
    end
end
