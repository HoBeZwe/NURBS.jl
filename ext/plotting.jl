
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
function NURBS.plotPatches(Patches; plotControlPoints=true, enforceRatio=true, resolution=0.05)

    uEvalpoints = collect(0:resolution:1.0)
    vEvalpoints = collect(0:resolution:1.0)

    traces = PlotlyJS.GenericTrace[]
    maxvec = Float64[]

    col = 1:size(Patches, 1)

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
                S; returnTrace=true, surfaceColor=col[ind] * ones(length(uEvalpoints), length(vEvalpoints)), cmin=1, cmax=maximum(col)
            )
        end
        [push!(traces, trace[i]) for i in eachindex(trace)]
        push!(maxvec, maxi)
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
