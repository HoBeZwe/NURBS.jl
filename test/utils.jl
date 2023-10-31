
@testset "Ranges" begin

    b = 9
    p = 3
    kVec = generateKnotVec(b, p)

    # example with points in every span
    evalpoints = 0:0.1:1
    Ranges = spanRanges(Bspline(p, kVec), evalpoints; emptyRanges=true)

    R = [0:-1, 0:-1, 0:-1, 1:2, 3:4, 5:5, 6:7, 8:9, 10:11, 0:-1, 0:-1, 0:-1]

    @test Ranges == R

    # example with spans without points
    evalpoints = [0.52]
    Ranges = spanRanges(Bspline(p, kVec), evalpoints; emptyRanges=true)

    R = [0:-1, 0:-1, 0:-1, 0:-1, 0:-1, 0:-1, 1:1, 0:-1, 0:-1, 0:-1, 0:-1, 0:-1]

    @test Ranges == R

    # example without empty ranges
    evalpoints = [0.52]
    Ranges = spanRanges(Bspline(p, kVec), evalpoints; emptyRanges=false)

    R = [0:-1, 0:-1, 0:-1, 1:1, 0:-1, 0:-1]

    @test Ranges == R
end

@testset "Greville + Anchors" begin

    b = 7
    p = 2
    kVec = generateKnotVec(b, p)

    Bspl = Bspline(p, kVec)

    gs = greville(Bspl)
    ac = anchors(Bspl)

    @test gs ≈ [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0] # computed by hand accrding to definitions
    @test ac ≈ [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0] # computed by hand accrding to definitions
end

@testset "File IO" begin

    @test_nowarn Patches = readMultipatch("assets/sphere.dat")

    Patches = readMultipatch("assets/sphere.dat")
    @test length(Patches) == 6
end

@testset "Plotting" begin

    @testset "plotPatches" begin

        Patches = readMultipatch("assets/sphere.dat")
        @test_nowarn plotPatches(Patches, plotControlPoints=true, resolution=0.1)
        @test_nowarn plotPatches(Patches, plotControlPoints=false, resolution=0.1, enforceRatio=false) # different optional arguments
    end

    @testset "Plotting Curves" begin

        kVec = Float64[0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 6, 6, 6]
        kVec ./= maximum(kVec)

        evalpoints = collect(0:0.05:1.0)

        p = 3               # degree of NURBS

        P1 = SVector(0.0, 0.0, 0.0)
        P2 = SVector(0.1, 0.25, 0.0)
        P3 = SVector(0.25, 0.3, 0.0)
        P4 = SVector(0.3, 0.5, 0.0)
        P5 = SVector(0.4, 0.4, 0.0)
        P6 = SVector(0.6, 0.3, 0.0)
        P7 = SVector(0.8, 0.7, 1.0)
        P8 = SVector(1.0, 0.4, 0.0)
        P9 = SVector(1.1, 0.4, 0.0)

        controlPoints = [P1, P2, P3, P4, P5, P6, P7, P8, P9]
        w             = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

        N = NURBScurve(NURB(p, kVec, w), controlPoints)
        C = N(evalpoints)

        @test_nowarn plotCurve3D(C, controlPoints=controlPoints)
        @test_nowarn plotCurve(C, controlPoints=controlPoints)

        C = N(evalpoints, 1)
        @test_nowarn plotCurve3D(C[1], controlPoints=controlPoints, tangents=C[2])
    end

    @testset "Plotting Surfaces" begin

        # --- parameters
        kVec = Float64[0, 0, 0, 1, 2, 3, 4, 5, 5, 5]
        kVec ./= maximum(kVec)

        p = 2               # degree of NURBS

        uEvalpoints = collect(0:0.01:1.0)
        vEvalpoints = collect(0:0.01:1.0)

        # ---
        controlPoints = [
            [
                SVector(0.0, 0.0, 3.0),
                SVector(0.0, 1.0, 3.0),
                SVector(0.0, 2.0, 2.0),
                SVector(0.0, 3.0, 2.0),
                SVector(0.0, 4.0, 2.0),
                SVector(0.0, 5.0, 2.0),
                SVector(0.0, 6.0, 2.0),
            ],
            [
                SVector(1.0, 0.0, 3.0),
                SVector(1.0, 1.0, 3.0),
                SVector(1.0, 2.0, 2.0),
                SVector(1.0, 3.0, 2.0),
                SVector(1.0, 4.0, 2.0),
                SVector(1.0, 5.0, 2.0),
                SVector(1.0, 6.0, 2.0),
            ],
            [
                SVector(2.0, 0.0, 2.0),
                SVector(2.0, 1.0, 2.0),
                SVector(2.0, 2.0, 1.0),
                SVector(2.0, 3.0, 1.0),
                SVector(2.0, 4.0, 1.0),
                SVector(2.0, 5.0, 1.0),
                SVector(2.0, 6.0, 1.0),
            ],
            [
                SVector(3.0, 0.0, 2.0),
                SVector(3.0, 1.0, 2.0),
                SVector(3.0, 2.0, 1.0),
                SVector(3.0, 3.0, 1.0),
                SVector(3.0, 4.0, 1.0),
                SVector(3.0, 5.0, 0.0),
                SVector(3.0, 7.0, 0.0),
            ],
            [
                SVector(4.0, 0.0, 1.0),
                SVector(4.0, 1.0, 1.0),
                SVector(4.0, 2.0, 0.0),
                SVector(4.0, 3.0, 0.0),
                SVector(4.0, 4.0, 1.0),
                SVector(4.0, 5.0, 0.0),
                SVector(4.0, 6.0, 0.0),
            ],
            [
                SVector(5.0, 0.0, 1.0),
                SVector(5.0, 1.0, 1.0),
                SVector(5.0, 2.0, 0.0),
                SVector(5.0, 3.0, 0.0),
                SVector(5.0, 4.0, 0.0),
                SVector(5.0, 5.0, 0.0),
                SVector(5.0, 6.0, 0.0),
            ],
            [
                SVector(6.0, 0.0, 1.0),
                SVector(6.0, 1.0, 1.0),
                SVector(6.0, 2.0, 0.0),
                SVector(6.0, 3.0, 0.0),
                SVector(6.0, 4.0, 0.0),
                SVector(6.0, 5.0, 0.0),
                SVector(6.0, 6.0, 0.0),
            ],
        ]

        controlPoints = [controlPoints[i][j] for i in 1:7, j in 1:7]

        Patch = BsplineSurface(Bspline(p, kVec), Bspline(p, kVec), controlPoints)

        # plot derivatives
        S = Patch(uEvalpoints, vEvalpoints, 1)

        @test_nowarn plotSurface(S[1, 1], tangents=S[2, 1])
        @test_nowarn plotSurface(S[1, 1], tangents=S[2, 1], controlPoints=Patch.controlPoints, enforceRatio=false)

        col = [1.0 for ui in eachindex(uEvalpoints), vi in eachindex(vEvalpoints)]
        @test_nowarn plotSurface(S[1, 1], surfaceColor=col)
    end
end
