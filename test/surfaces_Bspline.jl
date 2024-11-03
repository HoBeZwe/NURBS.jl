
using StaticArrays
using LinearAlgebra

@testset "B-splines" begin

    @testset "Surfaces Evaluation" begin

        # ------------------------------------------------ Example 3.4, p. 103 'The NURBS Book'
        # --- parameters
        p = 2 # degree 
        q = 2 # degree

        U = Float64[0, 0, 0, 2 / 5, 3 / 5, 1, 1, 1]
        V = Float64[0, 0, 0, 1 / 5, 1 / 2, 4 / 5, 1, 1, 1]

        uEvalpoints = [1 / 5]
        vEvalpoints = [3 / 5]

        # --- B-spline structure
        Bu = Bspline(p, U)
        Bv = Bspline(q, V)

        # --- resulting number of basis functions
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
        ]

        Ps = [controlPoints[i][j] for i in 1:5, j in 1:6]


        Patch = BsplineSurface(Bu, Bv, Ps)
        S = Patch(uEvalpoints, vEvalpoints)

        a = evalNaive(Bu, 1, uEvalpoints)[1]
        b = evalNaive(Bu, 2, uEvalpoints)[1]
        c = evalNaive(Bu, 3, uEvalpoints)[1]

        d = evalNaive(Bv, 3, vEvalpoints)[1]
        e = evalNaive(Bv, 4, vEvalpoints)[1]
        f = evalNaive(Bv, 5, vEvalpoints)[1]

        # actual analytical value
        Sref =
            d * (a * Ps[1, 3] + b * Ps[2, 3] + c * Ps[3, 3]) +
            e * (a * Ps[1, 4] + b * Ps[2, 4] + c * Ps[3, 4]) +
            f * (a * Ps[1, 5] + b * Ps[2, 5] + c * Ps[3, 5])

        @test Sref ≈ S[1, 1]
    end


    @testset "Derivatives Evaluation" begin

        # ------------------------------------------------ Example 3.4, p. 103 'The NURBS Book'
        # --- parameters
        p = 2 # degree 
        q = 2 # degree

        U = Float64[0, 0, 0, 2 / 5, 3 / 5, 1, 1, 1]
        V = Float64[0, 0, 0, 1 / 5, 1 / 2, 4 / 5, 1, 1, 1]

        uEvalpoints = [1 / 5]
        vEvalpoints = [3 / 5]

        # --- B-spline structure
        Bu = Bspline(p, U)
        Bv = Bspline(q, V)

        # --- resulting number of basis functions
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
        ]

        Ps = [controlPoints[i][j] for i in 1:5, j in 1:6]


        Patch = BsplineSurface(Bu, Bv, Ps)
        S = Patch(uEvalpoints, vEvalpoints, 1)

        a = evalNaive(Bu, 1, uEvalpoints)[1]
        b = evalNaive(Bu, 2, uEvalpoints)[1]
        c = evalNaive(Bu, 3, uEvalpoints)[1]

        d = evalNaive(Bv, 3, vEvalpoints)[1]
        e = evalNaive(Bv, 4, vEvalpoints)[1]
        f = evalNaive(Bv, 5, vEvalpoints)[1]

        # actual analytical value
        Sref =
            d * (a * Ps[1, 3] + b * Ps[2, 3] + c * Ps[3, 3]) +
            e * (a * Ps[1, 4] + b * Ps[2, 4] + c * Ps[3, 4]) +
            f * (a * Ps[1, 5] + b * Ps[2, 5] + c * Ps[3, 5])

        @test Sref ≈ S[1, 1][1, 1] # 0-th derivatives

        # the following three values are assumed to be correct in the current implementation
        @test S[1, 2][1, 1] ≈ SVector(0.0, 3.555555555, 0.0)
        @test S[2, 1][1, 1] ≈ SVector(4.1666666666, 0.0, -1.6666666666)
        @test isapprox(S[2, 2][1, 1], SVector(0.0, 0.0, 0.0), atol=1e-6)
    end
end
