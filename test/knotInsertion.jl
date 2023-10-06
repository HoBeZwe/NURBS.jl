
@testset "Knot insertion" begin

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
            kVecNew, ctrlPointsNew, wNew = insertKnot(kVec, controlPoints, p, uNew, 1)

            # evaluate curve with inserted point
            N = BsplineCurve(Bspline(p, kVecNew), ctrlPointsNew)
            C2 = N(evalpoints)

            # verify
            @test minimum(C1 .≈ C2)
        end

        @testset "Existing point twice" begin

            # insert already existing parametric point twice
            uNew = 3 / 5
            kVecNew, ctrlPointsNew, wNew = insertKnot(kVec, controlPoints, p, uNew, 2)

            # evaluate curve with inserted point
            N = BsplineCurve(Bspline(p, kVecNew), ctrlPointsNew)
            C2 = N(evalpoints)

            # verify
            @test minimum(C1 .≈ C2)
        end

        @testset "New point once" begin

            # insert already existing parametric point twice
            uNew = 3.5 / 5
            kVecNew, ctrlPointsNew, wNew = insertKnot(kVec, controlPoints, p, uNew, 1)

            # evaluate curve with inserted point
            N = BsplineCurve(Bspline(p, kVecNew), ctrlPointsNew)
            C2 = N(evalpoints)

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
            kVecNew, ctrlPointsNew, wNew = insertKnot(kVec, controlPoints, p, uNew, 1, w)

            # evaluate curve with inserted point
            N = NURBScurve(NURB(p, kVecNew, wNew), ctrlPointsNew)
            C2 = N(evalpoints)

            # verify
            @test minimum(C1 .≈ C2)
        end

        @testset "Existing point twice" begin

            # insert parametric point
            uNew = 3 / 5
            kVecNew, ctrlPointsNew, wNew = insertKnot(kVec, controlPoints, p, uNew, 2, w)

            # evaluate curve with inserted point
            N = NURBScurve(NURB(p, kVecNew, wNew), ctrlPointsNew)
            C2 = N(evalpoints)

            # verify
            @test minimum(C1 .≈ C2)
        end

        @testset "New point once" begin

            # insert parametric point
            uNew = 3.5 / 5
            kVecNew, ctrlPointsNew, wNew = insertKnot(kVec, controlPoints, p, uNew, 1, w)

            # evaluate curve with inserted point
            N = NURBScurve(NURB(p, kVecNew, wNew), ctrlPointsNew)
            C2 = N(evalpoints)

            # verify
            @test minimum(C1 .≈ C2)
        end
    end

end
