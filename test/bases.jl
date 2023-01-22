
@testset "B-splines" begin

    @testset "Basis Evaluation" begin

        # ------------------------------------------------ B-spline basis: partition of unity
        # --- parameters
        b = 6       # number of basis functions
        p = 2       # degree 

        evalpoints = collect(0:0.001:1.0)

        # --- resulting knot vector
        kVec = generateKnotVec(b, p)

        # --- B-spline structure
        Bspl = Bspline(p, kVec)

        # --- smart eval 
        Bsmart = bSpline(Bspl, evalpoints)

        # --- B-spline basis sums to 1 in each point
        @test minimum(sum(Bsmart; dims=2) .≈ 1.0)



        # ------------------------------------------------ Example 2.3, p. 68 'The NURBS Book'
        # --- parameters
        p = 2 # degree 

        kVec = Float64[0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 6, 6, 6]
        #kVec ./= maximum(kVec)

        evalpoints = [2.5 / 6]


        # --- B-spline structure
        Bspl = Bspline(p, kVec)

        # --- resulting number of basis functions
        @test numBasisFunctions(Bspl) == 10


        # --- naive eval
        B = bSplineNaive(Bspl, 4, evalpoints)
        @test B[1] ≈ 1 / 8
        @test Bspl.knotVec[end] == 1.0 # automatic normalization

        # --- smart eval 
        Bsmart = bSpline(Bspl, evalpoints)
        @test Bsmart[1] ≈ 1 / 8
    end


    @testset "Derivative Evaluation" begin

        # ------------------------------------------------ Example 2.4, p. 71 'The NURBS Book'
        # --- parameters
        p = 2 # degree 

        kVec = Float64[0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5]
        #kVec ./= maximum(kVec)

        evalpoints = [2.5 / 5]


        # --- B-spline structure
        Bspl = Bspline(p, kVec)


        # --- naive eval
        B = bSplineNaiveDerivative(Bspl, 4, 1, evalpoints)
        @test B[1] ≈ -2.5
        @test Bspl.knotVec[end] == 1.0 # automatic normalization

        B = bSplineNaiveDerivative(Bspl, 4, 2, evalpoints)
        @test B[1] ≈ 25.0

        # --- smart eval 
        Bsmart = bSplineDerivatives(Bspl, 2, evalpoints)
        @test Bsmart[2] ≈ -2.5
        @test Bsmart[3] ≈ 25.0
    end
end


@testset "NURBS" begin

    @testset "Basis Evaluation" begin

        # ------------------------------------------------ B-spline basis: partition of unity
        # --- parameters
        b = 6       # number of basis functions
        p = 2       # degree 

        w = ones(b)
        w[4] = 1.8

        evalpoints = collect(0:0.001:1.0)

        # --- resulting knot vector
        kVec = generateKnotVec(b, p)

        # --- B-spline structure
        Nrb = NURB(p, kVec, w)

        # --- eval
        N = nurbsNaive(Nrb, 1, evalpoints) .* 0.0
        for i in 1:6
            N .+= nurbsNaive(Nrb, i, evalpoints)
        end

        # --- B-spline basis sums to 1 in each point
        @test minimum(N[1:(end - 1)] .≈ 1.0)



        # ------------------------------------------------ Example 2.3, p. 68 'The NURBS Book' (setting all weights to 1)
        # --- parameters
        p = 2 # degree 

        kVec = Float64[0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 6, 6, 6]
        #kVec ./= maximum(kVec)

        w = ones(10)

        evalpoints = [2.5 / 6]


        # --- B-spline structure
        Nrb = NURB(p, kVec, w)

        # --- resulting number of basis functions
        @test numBasisFunctions(Nrb) == 10


        # --- naive eval
        N = nurbsNaive(Nrb, 4, evalpoints)
        @test N[1] ≈ 1 / 8
        @test Nrb.knotVec[end] == 1.0 # automatic normalization
    end

    @testset "Derivative Evaluation" begin

        # ------------------------------------------------ Example 2.4, p. 71 'The NURBS Book' (setting all weights to 1)
        # --- parameters
        p = 2 # degree 

        kVec = Float64[0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5]
        #kVec ./= maximum(kVec)

        w = ones(8)

        evalpoints = [2.5 / 5]


        # --- B-spline structure
        N = NURB(p, kVec, w)


        # --- naive eval
        B = nurbsNaiveDerivative(N, 4, 1, evalpoints)
        @test B[1] ≈ -2.5
        @test N.knotVec[end] == 1.0 # automatic normalization

        B = nurbsNaiveDerivative(N, 4, 2, evalpoints)
        @test B[1] ≈ 25.0
    end

end
