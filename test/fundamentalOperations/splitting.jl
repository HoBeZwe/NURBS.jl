
using StaticArrays
using LinearAlgebra

@testset "Splitting" begin

    @testset "Curves" begin

        kVec = Float64[0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5, 5]
        p = 4

        P1 = SVector(0.0, 0.0, 0.0)
        P2 = SVector(0.1, 0.25, 0.0)
        P3 = SVector(0.25, 0.3, 0.0)
        P4 = SVector(0.3, 0.5, 0.0)
        P5 = SVector(0.4, 0.4, 0.0)
        P6 = SVector(0.6, 0.3, 0.0)
        P7 = SVector(0.8, 0.7, 1.0)
        P8 = SVector(1.0, 0.4, 0.0)
        P8 = SVector(1.5, 0.4, 0.0)
        P9 = SVector(2.5, 0.4, 0.0)

        controlPoints = [P1, P2, P3, P4, P5, P6, P7, P8, P9]
        w             = [1.4, 0.1, 1.2, 1.5, 1.3, 1.9, 3.0, 2.0, 1.1]


        @testset "B-spline Curve" begin

            # --- original curve
            BC = BsplineCurve(Bspline(p, kVec), controlPoints)

            evPs = collect(0:0.002:1.0) # 501 points -> dividable by 3: 501/3 = 167
            Crv  = BC(evPs)

            # --- split the curve into 3 parts
            spCvec = split(BC, 3)
            evPs   = collect(0:0.006:1.0) # 167 points

            # --- test the first part
            @test Crv[1:167] ≈ spCvec[1](evPs)
            @test Crv[end] ≈ spCvec[3]([1.0])[end] # endpoint
        end

        @testset "NURBS Curve" begin

            # --- original curve
            NC = NURBScurve(NURB(p, kVec, w), controlPoints)

            @testset "Single Split" begin

                # --- eval original curve
                Crv = NC([0.0, 1.0])

                # --- split the curve
                spC1, spC2 = split(NC, 0.2)

                Crv1 = spC1([0.0, 1.0])
                Crv2 = spC2([0.0, 1.0])

                # --- test whether endpoints agree
                @test Crv[1] ≈ Crv1[1]
                @test Crv1[end] ≈ Crv2[1]
                @test Crv[end] ≈ Crv2[end]
            end

            @testset "Multiple Splits" begin

                # --- eval original curve
                evPs = collect(0:0.002:1.0) # 501 points -> dividable by 3: 501/3 = 167
                Crv  = NC(evPs)

                # --- split the curve into 3 parts
                spCvec = split(NC, 3)
                evPs   = collect(0:0.006:1.0) # 167 points

                # --- test the first part
                @test Crv[1:167] ≈ spCvec[1](evPs)
                @test Crv[end] ≈ spCvec[3]([1.0])[end] # endpoint
            end
        end
    end


    @testset "Surfaces" begin

        kVecU = Float64[0, 0, 0, 1, 2, 3, 4, 5, 5, 5]
        kVecV = Float64[0, 0, 0, 1, 2, 3, 5, 5, 5]

        p = 2               # degree of NURBS

        BsplU = Bspline(p, kVecU)
        BsplV = Bspline(p, kVecV)

        uEvalp = collect(0:0.01:1.0)
        vEvalp = collect(0:0.01:1.0)

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

        @testset "B-spline Surface" begin

            # --- original surface
            Patch = BsplineSurface(Bspline(p, kVecU), Bspline(p, kVecV), controlPoints)
            Srfc  = Patch(uEvalp, vEvalp)

            # --- split the surface into 9 parts
            sPvec = split(Patch, 3, 3)

            # --- test whether endpoints agree
            @test sPvec[1]([0.0], [0.0]) ≈ Patch([0.0], [0.0])
            @test sPvec[end]([1.0], [1.0]) ≈ Patch([1.0], [1.0])

            @test sPvec[1]([1.0], [0.0]) ≈ sPvec[4]([0.0], [0.0])
            @test sPvec[1]([0.0], [1.0]) ≈ sPvec[2]([0.0], [0.0])
            @test sPvec[1]([1.0], [1.0]) ≈ sPvec[5]([0.0], [0.0])
            @test sPvec[5]([1.0], [1.0]) ≈ sPvec[9]([0.0], [0.0])
        end

        @testset "NURBS Surface" begin

            # --- original surface
            Patch = NURBSsurface(Bspline(p, kVecU), Bspline(p, kVecV), controlPoints, w)
            Srfc  = Patch(uEvalp, vEvalp)

            # --- split the surface into 9 parts
            sPvec = split(Patch, 3, 3)

            # --- test whether endpoints agree
            @test sPvec[1]([0.0], [0.0]) ≈ Patch([0.0], [0.0])
            @test sPvec[end]([1.0], [1.0]) ≈ Patch([1.0], [1.0])

            @test sPvec[1]([1.0], [0.0]) ≈ sPvec[4]([0.0], [0.0])
            @test sPvec[1]([0.0], [1.0]) ≈ sPvec[2]([0.0], [0.0])
            @test sPvec[1]([1.0], [1.0]) ≈ sPvec[5]([0.0], [0.0])
            @test sPvec[5]([1.0], [1.0]) ≈ sPvec[9]([0.0], [0.0])
        end

    end
end
