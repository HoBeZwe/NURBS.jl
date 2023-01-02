
@testset "B-splines" begin

    @testset "Curves Evaluation" begin

        # ------------------------------------------------ Example 2.3, p. 81 'The NURBS Book'
        # --- parameters
        p = 2 # degree 

        kVec = Float64[0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5]
        kVec ./= maximum(kVec)

        evalpoints = [2.5 / 5]

        # --- B-spline structure
        Bspl = Bspline(p, kVec)

        # --- resulting number of basis functions
        P1 = SVector(0.0, 0.0, 0.0)
        P2 = SVector(0.1, 0.25, 0.0)
        P3 = SVector(0.25, 0.3, 0.0)
        P4 = SVector(0.3, 0.5, 5.0)
        P5 = SVector(0.4, 0.4, 0.0)
        P6 = SVector(0.6, 0.3, 0.0)
        P7 = SVector(0.8, 0.7, 1.0)
        P8 = SVector(1.0, 0.4, 0.0)


        controlPoints = [P1, P2, P3, P4, P5, P6, P7, P8]

        N = BsplineCurve(Bspl, controlPoints)
        C = curvePoints(N, evalpoints)

        # actual analytical value
        Cref = 1 / 8 * P3 + 6 / 8 * P4 + 1 / 8 * P5

        @test Cref ≈ C[1]
    end


    @testset "Derivatives Evaluation" begin

        # ------------------------------------------------ Example p. 91 'The NURBS Book'
        # --- parameters
        p = 2 # degree 

        kVec = Float64[0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5]
        kVec ./= maximum(kVec)

        evalpoints = [2.5 / 5]

        # --- B-spline structure
        Bspl = Bspline(p, kVec)

        # --- resulting number of basis functions
        P1 = SVector(0.0, 0.0, 0.0)
        P2 = SVector(0.1, 0.25, 0.0)
        P3 = SVector(0.25, 0.3, 0.0)
        P4 = SVector(0.3, 0.5, 5.0)
        P5 = SVector(0.4, 0.4, 0.0)
        P6 = SVector(0.6, 0.3, 0.0)
        P7 = SVector(0.8, 0.7, 1.0)
        P8 = SVector(1.0, 0.4, 0.0)


        controlPoints = [P1, P2, P3, P4, P5, P6, P7, P8]

        N = BsplineCurve(Bspl, controlPoints)
        C = curveDerivativesPoints(N, evalpoints, 1)

        # actual analytical value
        Cref = -5 / 2 * P3 + 5 / 2 * P5

        @test Cref ≈ C[2][1]
    end
end
