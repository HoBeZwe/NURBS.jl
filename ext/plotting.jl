
"""
    plotCurve3D(C; controlPoints=[], tangents=[], tangentRes=25)

Plot a 3D curve given the vector C containing an SVector(x, y, z) for each point of the curve.
"""
function NURBS.plotCurve3D(C; controlPoints=[], tangents=[], tangentRes=25, returnTrace=false)

    data = PlotlyJS.GenericTrace[]

    x = [C[i][1] for i in eachindex(C)]
    y = [C[i][2] for i in eachindex(C)]
    z = [C[i][3] for i in eachindex(C)]

    t1 = PlotlyJS.scatter3d(; x=x, y=y, z=z, mode="lines", markersize=0.4, legend=:none, showlegend=false)
    push!(data, t1)

    # plot tangent vectors?
    if !isempty(tangents)

        u = [tangents[i][1] for i in eachindex(tangents)]
        v = [tangents[i][2] for i in eachindex(tangents)]
        w = [tangents[i][3] for i in eachindex(tangents)]

        t3 = PlotlyJS.cone(;
            x=x[1:tangentRes:end],
            y=y[1:tangentRes:end],
            z=z[1:tangentRes:end],
            u=u[1:tangentRes:end],
            v=v[1:tangentRes:end],
            w=w[1:tangentRes:end],
            showscale=false,
            sizemode="scaled",
            sizeref=1.0,
        )
        push!(data, t3)
    end

    # plot the control points?
    if !isempty(controlPoints)
        maxmax, trace = plotControlPoints(controlPoints)
        push!(data, trace)
    end

    returnTrace && return data

    # plot all traces
    PlotlyJS.plot(data)
end


"""
    plotCurve(C; controlPoints=[])

Plot a 2D curve given the vector C containing an SVector(x, y, z) for each point of the curve.
"""
function NURBS.plotCurve(C; controlPoints=[])

    data = PlotlyJS.GenericTrace[]

    x = [C[i][1] for i in eachindex(C)]
    y = [C[i][2] for i in eachindex(C)]

    t1 = PlotlyJS.scatter(; x=x, y=y, mode="lines", markersize=0.4, legend=:none)
    push!(data, t1)

    # plot the control points?
    if !isempty(controlPoints)

        xP = [controlPoints[i][1] for i in eachindex(controlPoints)]
        yP = [controlPoints[i][2] for i in eachindex(controlPoints)]

        i = ["P$(i)" for i in eachindex(controlPoints)]

        t2 = PlotlyJS.scatter(;
            x=xP, y=yP, mode="lines+markers+text", markersize=0.4, legend=:none, text=i, textposition="top center", showlegend=false
        )

        push!(data, t2)
    end

    PlotlyJS.plot(data)
end


"""
    plotSurface(S; controlPoints=[], enforceRatio=true, returnTrace=false, surfaceColor=[], tangents=[])

Plot a 3D surface given the matrix 'S' containing an SVector(x, y, z) for each point of the surface.
"""
function NURBS.plotSurface(S; controlPoints=[], enforceRatio=true, returnTrace=false, surfaceColor=[], tangents=[], cmin=1, cmax=1)

    data = PlotlyJS.GenericTrace[]

    x = [S[i, j][1] for i in eachindex(S[:, 1]), j in eachindex(S[1, :])]
    y = [S[i, j][2] for i in eachindex(S[:, 1]), j in eachindex(S[1, :])]
    z = [S[i, j][3] for i in eachindex(S[:, 1]), j in eachindex(S[1, :])]

    if isempty(surfaceColor)
        col = z
        cmin = minimum(z)
        cmax = maximum(z)
    else
        col = surfaceColor
    end

    t1 = PlotlyJS.surface(;
        z=z, x=x, y=y, surfacecolor=col, intensitymode="cell", colorscale="Viridis", opacity=1.0, showscale=false, cmin=cmin, cmax=cmax
    )
    push!(data, t1)

    if !isempty(tangents)
        u = [tangents[i, j][1] for i in eachindex(tangents[:, 1]), j in eachindex(tangents[1, :])]
        v = [tangents[i, j][2] for i in eachindex(tangents[:, 1]), j in eachindex(tangents[1, :])]
        w = [tangents[i, j][3] for i in eachindex(tangents[:, 1]), j in eachindex(tangents[1, :])]

        t3 = PlotlyJS.cone(; x=x[:], y=y[:], z=z[:], u=u[:], v=v[:], w=w[:], showscale=false, sizemode="scaled", sizeref=1.0)
        push!(data, t3)
    end

    maxmax = 0.0

    # plot the control points?
    if !isempty(controlPoints)
        # along u-dirction
        for ind in eachindex(controlPoints[:, 1])
            maxmax, t2 = plotControlPoints(controlPoints[ind, :])
            push!(data, t2)
        end

        # along v-direction
        for ind in eachindex(controlPoints[1, :])
            maxmax2, t2 = plotControlPoints(controlPoints[:, ind])
            push!(data, t2)

            maxmax = maximum([maxmax, maxmax2])
        end
    end

    maxmax = maximum([maxmax, maximum(abs.(x)), maximum(abs.(y)), maximum(abs.(z))])

    if returnTrace
        return maxmax, data # do not plot but just return the traces
    end

    # ensure all three axes have the same length
    layout = PlotlyJS.Layout(;
        scene=PlotlyJS.attr(;
            xaxis=PlotlyJS.attr(; visible=true, legend=:none, range=[-maxmax, maxmax]),
            yaxis=PlotlyJS.attr(; visible=true, range=[-maxmax, maxmax]),
            zaxis=PlotlyJS.attr(; visible=true, range=[-maxmax, maxmax]),
            aspectratio=PlotlyJS.attr(; x=1, y=1, z=1),
        ),
    )

    if enforceRatio
        PlotlyJS.plot(data, layout)
    else
        PlotlyJS.plot(data)
    end
end


"""
    plotPatches(Patches; plotControlPoints=true, enforceRatio=true, resolution=0.05)

Plot multipatch geometry.
"""
function NURBS.plotPatches(
    Patches; mesh=[], pos=[], plotControlPoints=true, localVertices=false, patchID=false, enforceRatio=true, resolution=0.05, color=[]
)

    uEvalpoints = collect(0:resolution:1.0)
    vEvalpoints = collect(0:resolution:1.0)

    traces = PlotlyJS.GenericTrace[]
    maxvec = Float64[]

    col = 1 .+ rand(length(Patches))#1:size(Patches, 1)

    for (ind, patch) in enumerate(Patches)

        S = patch(uEvalpoints, vEvalpoints)

        if plotControlPoints
            maxi, trace = plotSurface(
                S;
                controlPoints=patch.controlPoints,
                returnTrace=true,
                surfaceColor=col[ind] * ones(length(uEvalpoints), length(vEvalpoints)),
                cmin=1,
                cmax=maximum(col),
            )
        else
            maxi, trace = plotSurface(
                S;
                returnTrace=true,
                surfaceColor=col[ind] * ones(length(uEvalpoints), length(vEvalpoints)),
                cmin=1,
                cmax=maximum(col),
                #S; returnTrace=true, surfaceColor=color[ind], cmin=0, cmax=2
            )
        end
        [push!(traces, trace[i]) for i in eachindex(trace)]
        push!(maxvec, maxi)

        isempty(mesh) || plotMesh!(traces, patch, mesh, uEvalpoints, vEvalpoints)
        localVertices && plotLocalvertices!(traces, patch)
        patchID && plotPatchID!(traces, patch, ind)
    end

    if !isempty(pos)
        x = [pos[i][1] for i in eachindex(pos)]
        y = [pos[i][2] for i in eachindex(pos)]
        z = [pos[i][3] for i in eachindex(pos)]
        txt = ["$i" for i in 1:length(pos)]
        push!(
            traces, PlotlyJS.scatter3d(; x=x, y=y, z=z, mode="markers+text", text=txt, markersize=0.4, legend=:none, showlegend=false)
        )
    end

    maxmax = maximum(maxvec)
    layout = PlotlyJS.Layout(;
        scene=PlotlyJS.attr(;
            xaxis=PlotlyJS.attr(; visible=true, legend=:none, range=[-maxmax, maxmax]),
            yaxis=PlotlyJS.attr(; visible=true, range=[-maxmax, maxmax]),
            zaxis=PlotlyJS.attr(; visible=true, range=[-maxmax, maxmax]),
            aspectratio=PlotlyJS.attr(; x=1, y=1, z=1),
        ),
    )

    if enforceRatio
        PlotlyJS.plot(traces, layout)
    else
        PlotlyJS.plot(traces)
    end
end

""" 
    plotPatchID!(traces, patch, ind)

Plot the unique IDs of the patches.
"""
function plotPatchID!(traces, patch, ind)

    p1 = patch([0.5], [0.5])[1]

    t1 = PlotlyJS.scatter3d(;
        x=[p1[1]], y=[p1[2]], z=[p1[3]], mode="markers+text", text=["$ind"], markersize=0.4, legend=:none, showlegend=false
    )
    push!(traces, t1)

    return nothing
end


"""
    plotLocalvertices!(traces, patch)

Plot the four local vertices on each patch.
"""
function plotLocalvertices!(traces, patch)

    off = 0.025
    p1 = patch([off], [off])
    p2 = patch([off], [1 - off])
    p3 = patch([1 - off], [1 - off])
    p4 = patch([1 - off], [off])

    x = [p1[1][1], p2[1][1], p3[1][1], p4[1][1]]
    y = [p1[1][2], p2[1][2], p3[1][2], p4[1][2]]
    z = [p1[1][3], p2[1][3], p3[1][3], p4[1][3]]

    t1 = PlotlyJS.scatter3d(;
        x=x, y=y, z=z, mode="markers+text", text=["1", "2", "3", "4"], markersize=0.7, legend=:none, showlegend=false
    )
    push!(traces, t1)

    return nothing
end


"""
    plotMesh!(traces, patch, mesh, uEvalpoints, vEvalpoints)

Plot the mesh structure defined by mesh=[U, V] with U and V being vectors with entries in [0,1].
"""
function plotMesh!(traces, patch, mesh, uEvalpoints, vEvalpoints)

    uKnotMesh = mesh[1]
    vKnotMesh = mesh[2]

    for (ind, vVal) in enumerate(vKnotMesh)

        C = patch(uEvalpoints, [vVal])

        x = [C[i][1] for i in eachindex(C)]
        y = [C[i][2] for i in eachindex(C)]
        z = [C[i][3] for i in eachindex(C)]

        t1 = PlotlyJS.scatter3d(; x=x, y=y, z=z, mode="lines", markersize=0.4, legend=:none, showlegend=false)
        push!(traces, t1)
    end

    for (ind, uVal) in enumerate(uKnotMesh)

        C = patch([uVal], vEvalpoints)

        x = [C[i][1] for i in eachindex(C)]
        y = [C[i][2] for i in eachindex(C)]
        z = [C[i][3] for i in eachindex(C)]

        t1 = PlotlyJS.scatter3d(; x=x, y=y, z=z, mode="lines", markersize=0.4, legend=:none, showlegend=false)
        push!(traces, t1)
    end

    return nothing
end


"""
    plotControlPoints(controlPoints)

Plot points in 3D space and connect them with lines. 
The plot is not directly generated but the trace is returned.
"""
function plotControlPoints(controlPoints)

    xP = [controlPoints[i][1] for i in eachindex(controlPoints)]
    yP = [controlPoints[i][2] for i in eachindex(controlPoints)]
    zP = [controlPoints[i][3] for i in eachindex(controlPoints)]

    i = ["P$(i)" for i in eachindex(controlPoints)]

    maxmax = maximum([maximum(abs.(xP)), maximum(abs.(yP)), maximum(abs.(zP))])

    return maxmax,
    PlotlyJS.scatter3d(;
        x=xP, y=yP, z=zP, mode="lines+markers", markersize=0.4, legend=:none, textposition="top center", showlegend=false
    )
end
