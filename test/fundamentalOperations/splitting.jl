
@testset "Curve Splitting" begin

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
