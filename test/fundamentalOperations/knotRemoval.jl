
using StaticArrays
using LinearAlgebra

@testset "Knot insertion" begin

    @testset "Curves" begin

        # --- example from Rogers "An Introduction to NURBS" p. 117
        p = 3

        P1 = SVector(0.0, 0.0, 0.0)
        P2 = SVector(0.0, 2.0, 0.0)
        P3 = SVector(1.5, 3.0, 0.0)
        P4 = SVector(3.0, 3.0, 0.0)
        P5 = SVector(4.5, 3.0, 0.0)
        P6 = SVector(6.0, 2.0, 0.0)
        P7 = SVector(6.0, 0.0, 0.0)

        kVec = Float64[0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2]

        cP = [P1, P2, P3, P4, P5, P6, P7]
        w  = [1.4, 0.1, 1.2, 1.5, 1.3, 1.9, 3.0]

        evalpoints = collect(0:0.005:1.0)

        @testset "B-splines" begin

            BC = BsplineCurve(Bspline(p, kVec), cP)

            RC = removeKnot(BC, 0.5, 3)
            kVecMod, cPmod, wMod = coarsen(kVec, cP, p, [0.5, 0.5, 0.5])

            @test RC.basis.knotVec == kVecMod
            @test RC.controlPoints == cPmod

            cPR = RC.controlPoints

            @test cPR[1] ≈ SVector(0.0, 0.0, 0.0)
            @test cPR[2] ≈ SVector(0.0, 4.0, 0.0)
            @test cPR[3] ≈ SVector(6.0, 4.0, 0.0)
            @test cPR[4] ≈ SVector(6.0, 0.0, 0.0)

            C1 = BC(evalpoints)
            C2 = RC(evalpoints)

            @test C1 ≈ C2

            # --- remove knot more often than existent
            @test_logs (:warn, "The knot is only contained 3 times.") removeKnot(BC, 0.5, 4)

            # --- change curve such that 0.5 is only 1 time removable
            BC.controlPoints[6] = SVector(6.0, 2.0, 3.0)
            @test_logs (:warn, "The knot is only 1 times removable.") removeKnot(BC, 0.5, 2)
        end

        @testset "NURBS" begin

            cP[5] = SVector(4.5, 3.0, 2.1)
            cP[3] = SVector(1.5, 3.0, -1.0) # 3D example

            NC = NURBScurve(NURB(p, kVec, w), cP)

            IC = insertKnot(NC, 0.25, 2)
            IC = insertKnot(IC, 0.75, 3)

            RC = removeKnot(IC, 0.25, 1) # remove some of the inserted knots
            RC = removeKnot(RC, 0.75, 2) # remove some of the inserted knots

            C1 = NC(evalpoints)
            C2 = RC(evalpoints)

            @test C1 ≈ C2


            # --- try to remove endpoint knots
            @test_throws ErrorException("The knot has to be an interior knot.") removeKnot(NC, 0.0, 1)

            # --- remove non-existent knot
            @test_logs (:warn, "The requested knot is not contained in the knot vector.") removeKnot(NC, 0.33, 1)
        end
    end

    @testset "Surfaces" begin

        # --- surface with 3 times removable knot 0.5
        p = 3
        kVec = Float64[0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2]

        P1 = SVector(0.0, 0.0, 0.0)
        P2 = SVector(0.0, 2.0, 0.0)
        P3 = SVector(1.5, 3.0, 0.0)
        P4 = SVector(3.0, 3.0, 0.0)
        P5 = SVector(4.5, 3.0, 0.0)
        P6 = SVector(6.0, 2.0, 0.0)
        P7 = SVector(6.0, 0.0, 0.0)

        cP = [P1, P2, P3, P4, P5, P6, P7]

        P1 = SVector(0.0, 0.0, 2.0)
        P2 = SVector(0.0, 2.0, 2.0)
        P3 = SVector(1.5, 3.0, 2.0)
        P4 = SVector(3.0, 3.0, 2.0)
        P5 = SVector(4.5, 3.0, 2.0)
        P6 = SVector(6.0, 2.0, 2.0)
        P7 = SVector(6.0, 0.0, 2.0)

        cP2 = [P1, P2, P3, P4, P5, P6, P7]

        uEval = vEval = collect(0:0.005:1.0)

        @testset "B-splines" begin

            # --- in u-direction
            controlPoints = [[cP, cP2][j][i] for i in 1:7, j in 1:2]

            BS1 = BsplineSurface(Bspline(p, kVec), Bspline(1, [0.0, 0.0, 1.0, 1.0]), controlPoints)

            BS2 = removeKnotU(BS1, 0.5, 1) # remove knot once

            S1 = BS1(uEval, vEval)
            S2 = BS2(uEval, vEval)

            @test S1 ≈ S2 # before and after removal identical        


            # --- try to remove more often than possible
            BS1.controlPoints[4, 1] = SVector(5.5, 0.0, 6.6)
            @test_logs (:warn, "The knot 0.5 can not be removed 4 times.") match_mode = :any removeKnotU(BS1, 0.5, 4)


            # --- in v-direction
            controlPoints = [[cP, cP2][i][j] for i in 1:2, j in 1:7] # along v

            BS1 = BsplineSurface(Bspline(1, [0.0, 0.0, 1.0, 1.0]), Bspline(p, kVec), controlPoints)

            BS2 = removeKnotV(BS1, 0.5, 1)

            S1 = BS1(uEval, vEval)
            S2 = BS2(uEval, vEval)

            @test S1 ≈ S2 # before and after removal identical 


            # --- try to remove more often than possible
            BS1.controlPoints[1, 4] = SVector(5.5, 0.0, 6.6)
            @test_logs (:warn, "The knot 0.5 can not be removed 4 times.") match_mode = :any removeKnotV(BS1, 0.5, 4)
        end

        @testset "NURBS" begin

            controlPoints = [[cP, cP2][j][i] for i in 1:7, j in 1:2]

            w = [1.4, 0.1, 1.2, 1.5, 1.3, 1.9, 3.0]
            w = [[w, w][j][i] for i in 1:7, j in 1:2]

            NS1 = NURBSsurface(Bspline(p, kVec), Bspline(1, [0.0, 0.0, 1.0, 1.0]), controlPoints, w)

            NS2 = refine(NS1; U=[0.6, 0.6]) # insert a point twice

            # --- remove knot
            NS3 = removeKnotU(NS2, 0.6, 1) # remove the knot once

            S1 = NS2(uEval, vEval)
            S2 = NS3(uEval, vEval)

            @test S1 ≈ S2
        end
    end
end
