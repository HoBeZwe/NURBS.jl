
using StaticArrays
using LinearAlgebra

@testset "Knot insertion" begin

    @testset "Curves" begin

        # define setting
        kVec = Float64[0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5]
        kVec ./= maximum(kVec)

        p = 3 # polynomial degree

        P1 = SVector(0.0, 0.0, 0.0)
        P2 = SVector(0.1, 0.25, 0.0)
        P3 = SVector(0.25, 0.3, 0.0)
        P4 = SVector(0.3, 0.5, 0.0)
        P5 = SVector(0.4, 0.4, 0.0)
        P6 = SVector(0.6, 0.3, 0.0)
        P7 = SVector(0.8, 0.7, 1.0)
        P8 = SVector(1.0, 0.4, 0.0)
        P8 = SVector(1.5, 0.4, 0.0)

        controlPoints = [P1, P2, P3, P4, P5, P6, P7, P8]
        w             = [1.0, 0.1, 1.0, 1.0, 1.0, 1.0, 3.0, 1.0]

        evalpoints = collect(0:0.05:1.0)


        @testset "B-splines" begin

            # evaluate original curve
            N  = BsplineCurve(Bspline(p, kVec), controlPoints)
            C1 = N(evalpoints)

            @testset "Existing point once" begin

                # insert already existing parametric point twice
                uNew = 3 / 5
                N2 = insertKnot(N, uNew, 1)

                C2 = N2(evalpoints)

                # verify
                @test minimum(C1 .≈ C2)
            end

            @testset "Existing point twice" begin

                # insert already existing parametric point twice
                uNew = 3 / 5
                N2 = insertKnot(N, uNew, 2)

                C2 = N2(evalpoints)

                # verify
                @test minimum(C1 .≈ C2)
            end

            @testset "New point once" begin

                # insert new parametric point once
                uNew = 3.5 / 5
                N2 = insertKnot(N, uNew, 1)

                C2 = N2(evalpoints)

                # verify
                @test minimum(C1 .≈ C2)
            end

            @testset "Refine" begin

                # insert a list of parametric points
                N2 = refine(N, [0.1, 0.2, 0.3, 0.8221])
                C2 = N2(evalpoints)

                # verify
                @test minimum(C1 .≈ C2)
            end

            @testset "Degree = 0 case" begin

                kV2 = Float64[0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
                cP = [P1, P2, P3, P4, P5]

                N  = BsplineCurve(Bspline(0, kV2), cP)
                C1 = N(evalpoints)

                # insert a list of parametric points
                N2 = insertKnot(N, 0.1, 2)
                C2 = N2(evalpoints)

                # verify
                @test minimum(C1 .≈ C2)
            end
        end


        @testset "NURBS" begin

            # evaluate original curve
            N  = NURBScurve(NURB(p, kVec, w), controlPoints)
            C1 = N(evalpoints)

            @testset "Existing point once" begin

                # insert parametric point
                uNew = 3 / 5
                N2 = insertKnot(N, uNew, 1)

                C2 = N2(evalpoints)

                # verify
                @test minimum(C1 .≈ C2)
            end

            @testset "Existing point twice" begin

                # insert parametric point
                uNew = 3 / 5
                N2 = insertKnot(N, uNew, 2)

                C2 = N2(evalpoints)

                # verify
                @test minimum(C1 .≈ C2)
            end

            @testset "New point once" begin

                # insert parametric point
                uNew = 3.5 / 5
                N2 = insertKnot(N, uNew, 1)

                C2 = N2(evalpoints)

                # verify
                @test minimum(C1 .≈ C2)
            end

            @testset "New point more often than allowed" begin

                # insert parametric point
                uNew = 4.2 / 5
                N2 = insertKnot(N, uNew, p + 1)

                C2 = N2(evalpoints)

                # verify
                @test minimum(C1 .≈ C2)
            end

            @testset "Refine" begin

                # insert parametric point
                N2 = refine(N, [0.1, 0.2, 0.3, 0.8221])
                C2 = N2(evalpoints)

                # verify
                @test minimum(C1 .≈ C2)
            end
        end
    end


    @testset "Surfaces" begin

        kVecU = Float64[0, 0, 0, 1, 2, 3, 4, 5, 5, 5]
        kVecV = Float64[0, 0, 0, 1, 2, 3, 5, 5, 5]

        p = 2               # degree of NURBS

        BsplU = Bspline(p, kVecU)
        BsplV = Bspline(p, kVecV)

        uEvalpoints = collect(0:0.01:1.0)
        vEvalpoints = collect(0:0.01:1.0)

        controlPoints = [
            [
                SVector(0.0, 0.0, 3.0),
                SVector(0.0, 1.0, 3.0),
                SVector(0.0, 2.0, 2.0),
                SVector(0.0, 3.0, 2.0),
                SVector(0.0, 4.0, 2.0),
                SVector(0.0, 5.0, 2.0),
            ],
            [
                SVector(1.0, 0.0, 3.0),
                SVector(1.0, 1.0, 3.0),
                SVector(1.0, 2.0, 2.0),
                SVector(1.0, 3.0, 2.0),
                SVector(1.0, 4.0, 2.0),
                SVector(1.0, 5.0, 2.0),
            ],
            [
                SVector(2.0, 0.0, 2.0),
                SVector(2.0, 1.0, 2.0),
                SVector(2.0, 2.0, 1.0),
                SVector(2.0, 3.0, 1.0),
                SVector(2.0, 4.0, 1.0),
                SVector(2.0, 5.0, 1.0),
            ],
            [
                SVector(3.0, 0.0, 2.0),
                SVector(3.0, 1.0, 2.0),
                SVector(3.0, 2.0, 1.0),
                SVector(3.0, 3.0, 1.0),
                SVector(3.0, 4.0, 1.0),
                SVector(3.0, 5.0, 0.0),
            ],
            [
                SVector(4.0, 0.0, 1.0),
                SVector(4.0, 1.0, 1.0),
                SVector(4.0, 2.0, 0.0),
                SVector(4.0, 3.0, 0.0),
                SVector(4.0, 4.0, 1.0),
                SVector(4.0, 5.0, 0.0),
            ],
            [
                SVector(5.0, 0.0, 1.0),
                SVector(5.0, 1.0, 1.0),
                SVector(5.0, 2.0, 0.0),
                SVector(5.0, 3.0, 0.0),
                SVector(5.0, 4.0, 0.0),
                SVector(5.0, 5.0, 0.0),
            ],
            [
                SVector(6.0, 0.0, 1.0),
                SVector(6.0, 1.0, 1.0),
                SVector(6.0, 2.0, 0.0),
                SVector(6.0, 3.0, 0.0),
                SVector(6.0, 4.0, 0.0),
                SVector(6.0, 5.0, 0.0),
            ],
        ]

        controlPoints = [controlPoints[i][j] for i in 1:7, j in 1:6]

        w = ones(size(controlPoints))
        w[5, 5] = 5.0
        w[7, 2] = 0.1
        w[2, 1] = 7.0

        @testset "B-splines" begin

            Patch = BsplineSurface(Bspline(p, kVecU), Bspline(p, kVecV), controlPoints)

            PatchRef = refine(Patch; U=[0.5, 0.5, 0.72], V=[0.2, 0.42])

            P1 = Patch(uEvalpoints, vEvalpoints)
            P2 = PatchRef(uEvalpoints, vEvalpoints)

            @test P1 ≈ P2
        end

        @testset "NURBS" begin

            Patch = NURBSsurface(Bspline(p, kVecU), Bspline(p, kVecV), controlPoints, w)

            PatchRef = refine(Patch; U=[0.5, 0.5, 0.72], V=[0.2, 0.42])

            P1 = Patch(uEvalpoints, vEvalpoints)
            P2 = PatchRef(uEvalpoints, vEvalpoints)

            @test P1 ≈ P2
        end
    end
end
