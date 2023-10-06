
@testset "NURBS" begin

    @testset "Curves Evaluation" begin

        # ------------------------------------------------ Example 4.1, p. 122 'The NURBS Book'
        # --- parameters
        p = 2 # degree 

        kVec = Float64[0, 0, 0, 1, 2, 3, 3, 3]
        kVec ./= maximum(kVec)

        evalpoints = [1 / 3]

        w = Float64[1, 4, 1, 1, 1]

        # --- B-spline structure
        Nrb = NURB(p, kVec, w)

        # --- resulting number of basis functions
        P1 = SVector(0.0, 0.0, 0.0)
        P2 = SVector(1.0, 1.0, 0.0)
        P3 = SVector(3.0, 2.0, 0.0)
        P4 = SVector(4.0, 1.0, 0.0)
        P5 = SVector(5.0, -1.0, 0.0)

        controlPoints = [P1, P2, P3, P4, P5]

        N = NURBScurve(Nrb, controlPoints)
        C = N(evalpoints)

        # actual analytical value
        Cref = SVector(7 / 5, 6 / 5, 0.0)

        @test Cref ≈ C[1]
    end


    @testset "Derivatives Evaluation" begin

        # ------------------------------------------------ Example 4.2, p. 126 'The NURBS Book'
        # --- parameters
        p = 2 # degree 

        kVec = Float64[0, 0, 0, 1, 1, 1]
        kVec ./= maximum(kVec)

        evalpoints = [0.0, 1.0]

        w = Float64[1, 1, 2]

        # --- B-spline structure
        Nrb = NURB(p, kVec, w)

        # --- resulting number of basis functions
        P1 = SVector(1.0, 0.0, 0.0)
        P2 = SVector(1.0, 1.0, 0.0)
        P3 = SVector(0.0, 1.0, 0.0)

        controlPoints = [P1, P2, P3]

        N = NURBScurve(Nrb, controlPoints)
        C = N(evalpoints, 1)

        # actual analytical value
        Cref1 = SVector(0.0, 2.0, 0.0)
        Cref2 = SVector(-1.0, 0.0, 0.0)

        @test Cref1 ≈ C[2][1]
        @test Cref2 ≈ C[2][2]
    end
end
