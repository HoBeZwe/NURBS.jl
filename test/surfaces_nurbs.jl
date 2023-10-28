
@testset "NURBS" begin

    @testset "Surfaces Evaluation" begin

        # ------------------------------------------------ Example 4.3, p. 132 'The NURBS Book'
        # --- parameters
        p = 2 # degree 
        q = 2 # degree

        U = Float64[0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5]
        V = Float64[0, 0, 0, 1, 2, 3, 3, 3]

        U ./= 5.0
        V ./= 3.0

        uEvalpoints = [2.5 / 5]
        vEvalpoints = [1 / 3]

        # --- B-spline structure
        Bu = Bspline(p, U)
        Bv = Bspline(q, V)

        numBasisFunctions(Bu)
        numBasisFunctions(Bv)

        # --- resulting number of basis functions
        controlPoints = [
            [SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0)],
            [SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0)],
            [SVector(0.0, 0.0, 0.0), SVector(0.0, 2.0, 4.0), SVector(0.0, 3.0, 2.0), SVector(0.0, 2.0, 0.0), SVector(0.0, 0.0, 0.0)],
            [SVector(0.0, 0.0, 0.0), SVector(2.0, 3.0, 4.0), SVector(2.0, 4.0, 2.0), SVector(2.0, 3.0, 0.0), SVector(0.0, 0.0, 0.0)],
            [SVector(0.0, 0.0, 0.0), SVector(4.0, 2.0, 4.0), SVector(4.0, 3.0, 2.0), SVector(4.0, 2.0, 0.0), SVector(0.0, 0.0, 0.0)],
            [SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0)],
            [SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0)],
            [SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0)],
        ]

        Ps = [controlPoints[i][j] for i in 1:8, j in 1:5]
        w = [
            0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0
            0.0 1.0 2.0 1.0 0.0
            0.0 2.0 6.0 2.0 0.0
            0.0 1.0 2.0 1.0 0.0
            0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0
        ]

        Patch = NURBSsurface(Bu, Bv, Ps, w)
        S = Patch(uEvalpoints, vEvalpoints)

        Sref = SVector(2.0, 98 / 27, 68 / 27)

        @test Sref ≈ S[1, 1]
    end


    @testset "Derivatives Evaluation" begin

        # ------------------------------------------------ Example 3.4, p. 103 'The NURBS Book' (setting all w=1)
        # --- parameters
        p = 2 # degree 
        q = 2 # degree

        U = Float64[0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5]
        V = Float64[0, 0, 0, 1, 2, 3, 3, 3]

        U ./= 5.0
        V ./= 3.0

        uEvalpoints = [2.5 / 5]
        vEvalpoints = [1 / 3]

        # --- B-spline structure
        Bu = Bspline(p, U)
        Bv = Bspline(q, V)

        numBasisFunctions(Bu)
        numBasisFunctions(Bv)

        # --- resulting number of basis functions
        controlPoints = [
            [SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0)],
            [SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0)],
            [SVector(0.0, 0.0, 0.0), SVector(0.0, 2.0, 4.0), SVector(0.0, 3.0, 2.0), SVector(0.0, 2.0, 0.0), SVector(0.0, 0.0, 0.0)],
            [SVector(0.0, 0.0, 0.0), SVector(2.0, 3.0, 4.0), SVector(2.0, 4.0, 2.0), SVector(2.0, 3.0, 0.0), SVector(0.0, 0.0, 0.0)],
            [SVector(0.0, 0.0, 0.0), SVector(4.0, 2.0, 4.0), SVector(4.0, 3.0, 2.0), SVector(4.0, 2.0, 0.0), SVector(0.0, 0.0, 0.0)],
            [SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0)],
            [SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0)],
            [SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0)],
        ]

        Ps = [controlPoints[i][j] for i in 1:8, j in 1:5]
        w = [
            0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0
            0.0 1.0 2.0 1.0 0.0
            0.0 2.0 6.0 2.0 0.0
            0.0 1.0 2.0 1.0 0.0
            0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0
        ]

        Patch = NURBSsurface(Bu, Bv, Ps, w)
        S = Patch(uEvalpoints, vEvalpoints, 1)

        Sref = SVector(2.0, 98 / 27, 68 / 27) # actual analytical value

        @test Sref ≈ S[1, 1][1, 1] # no derivative

        # the following three values are assumed to be correct in the current implementation
        @test S[1, 2][1, 1] ≈ SVector(0.0, 2.403292181, -4.6090534979)
        @test S[2, 1][1, 1] ≈ SVector(4.4444444444, 0.0, 0.0)
        @test S[2, 2][1, 1] ≈ SVector(-3.950617283, 0.0, 0.0)

        @test_nowarn Jacobian(Patch, uEvalpoints, vEvalpoints)




        # ---- with memory preallocation
        pM = NURBS.preAllocNURBSsurface(p, q, uEvalpoints, vEvalpoints, 1)

        S = Patch(uEvalpoints, vEvalpoints, 1, pM)

        @test Sref ≈ S[1, 1][1, 1] # no derivative

        # the following three values are assumed to be correct in the current implementation
        @test S[1, 2][1, 1] ≈ SVector(0.0, 2.403292181, -4.6090534979)
        @test S[2, 1][1, 1] ≈ SVector(4.4444444444, 0.0, 0.0)
        @test S[2, 2][1, 1] ≈ SVector(-3.950617283, 0.0, 0.0)
    end
end
