
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

@testset "Convenience Functions" begin

    # --- normalization of knot vectors
    kVec = Float64[-1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 6, 6, 6]
    NURBS.isValidKnotVector!(kVec)

    @test kVec[1] == 0
    @test kVec[end] == 1

    # --- invalid knot vector
    kVec = Float64[-1, -1, -1, -1, 0, 1, 4, 3, 2, 5, 6, 6, 6, 6]
    @test_throws ErrorException("The knot vector is not in ascending order.") NURBS.isValidKnotVector!(kVec)

    # --- return weights
    w = ones(9)
    kVec = Float64[0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 6, 6, 6]
    NURBS.weights(NURB(3, kVec, w)) == w
end

@testset "Geometry Manipulation" begin

    Patches = readMultipatch("assets/sphere.dat")

    @testset "Scaling" begin

        fac = 2.3
        P2 = scale(Patches, fac)

        norms = norm.(P2[1].controlPoints)

        # length of the edge points of the patches correspond to the new radius
        @test norms[1, 1] ≈ fac
        @test norms[1, end] ≈ fac
        @test norms[end, 1] ≈ fac
        @test norms[end, end] ≈ fac
    end

    @testset "Translation" begin

        shift = SVector(2.3, 1.2, -4.2)
        P2 = translate(Patches, shift)

        for (i, p2) in enumerate(P2)
            @test p2.controlPoints[1, 1] ≈ Patches[i].controlPoints[1, 1] + shift
        end
    end

    @testset "Rotation" begin

        rotAxis = SVector(0.0, 0.0, 1.0)
        α = π / 2

        P1 = rotate(Patches[1], rotAxis, α) # 90° rotation around z-axis
        @test Patches[2].controlPoints ≈ P1.controlPoints # second patch agrees with the rotated one

        P1 = rotate(Patches, rotAxis, 2π) # 360° rotation around z-axis
        @test Patches[2].controlPoints ≈ P1[2].controlPoints # second patch agrees with the rotated one
    end

    @testset "Mirroring" begin

        normal = SVector(0.0, 1.0, 0.0)
        anchor = SVector(0.0, 0.0, 0.0)

        P2 = mirror(Patches, normal, anchor) # mirror in xz-plane

        @test P2[3].controlPoints[1, end] ≈ Patches[1].controlPoints[1, end]
    end
end

@testset "Greville + Anchors" begin
    b = 7
    p = 2
    kVec = generateKnotVec(b, p)

    Bspl = Bspline(p, kVec)

    for spline in [Bspline(p, kVec), NURB(p, kVec, ones(length(kVec) - p - 1)), CurrySchoenberg(p, kVec)]
        gs = greville(spline)
        ac = anchors(spline)

        @test gs ≈ [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0] # computed by hand according to definitions
        @test ac ≈ [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0] # computed by hand according to definitions
    end
end

@testset "File IO" begin

    # --- multipatch files
    @test_nowarn Patches = readMultipatch("assets/sphere.dat")

    Patches = readMultipatch("assets/sphere.dat")
    @test length(Patches) == 6

    Patches = load("assets/sphere.dat")
    @test length(Patches) == 6


    # --- .step files
    # B-spline
    @test_nowarn Patches = readStep("assets/torus.stp")

    Patches = readStep("assets/torus.stp")
    @test length(Patches) == 30

    Patches = load("assets/torus.stp")
    @test length(Patches) == 30

    # NURBS
    @test_nowarn Patches = readStep("assets/sphere.stp")

    Patches = readStep("assets/sphere.stp")
    @test length(Patches) == 6


    # --- .vtk files
    @test_nowarn Patches = saveVtk("test", Patches; resolution=0.02)
end

@testset "Plotting" begin

    @testset "plotPatches" begin

        Patches = readMultipatch("assets/sphere.dat")
        @test_nowarn plotPatches(Patches, plotControlPoints=true, resolution=0.25)
        @test_nowarn plotPatches(Patches, plotControlPoints=false, resolution=0.25, enforceRatio=false) # different optional arguments
        @test_nowarn plotPatches(Patches, plotControlPoints=false, resolution=0.25, localVertices=true, patchID=true) # different optional arguments
        @test_nowarn plotPatches(
            Patches, plotControlPoints=false, resolution=0.25, mesh=[[0.2, 0.1], [0.1, 0.7]], pos=[SVector(0.1, 0.5, 0.9)]
        )
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
