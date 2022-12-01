
export plotCurve3D, plotCurve

function __init__()
    @require PlotlyJS = "f0f68f2c-4968-5e81-91da-67840de0976a" begin

        @eval function plotCurve3D(C, controlPoints; ctrlPoints=true)

            data = PlotlyJS.GenericTrace[]

            x = [C[i][1] for i in eachindex(C)]
            y = [C[i][2] for i in eachindex(C)]
            z = [C[i][3] for i in eachindex(C)]

            t1 = PlotlyJS.scatter3d(; x=x, y=y, z=z, mode="lines", markersize=0.4, legend=:none)

            push!(data, t1)

            # plot the control points?
            if ctrlPoints

                xP = [controlPoints[i][1] for i in eachindex(controlPoints)]
                yP = [controlPoints[i][2] for i in eachindex(controlPoints)]
                zP = [controlPoints[i][3] for i in eachindex(controlPoints)]

                i = ["P$(i-1)" for i in eachindex(controlPoints)]

                t2 = PlotlyJS.scatter3d(;
                    x=xP, y=yP, z=zP, mode="lines+markers+text", markersize=0.4, legend=:none, text=i, textposition="top center"
                )

                push!(data, t2)
            end

            PlotlyJS.plot(data)
        end

        @eval function plotCurve(C, controlPoints; ctrlPoints=true)

            data = PlotlyJS.GenericTrace[]

            x = [C[i][1] for i in eachindex(C)]
            y = [C[i][2] for i in eachindex(C)]

            t1 = PlotlyJS.scatter(; x=x, y=y, mode="lines", markersize=0.4, legend=:none)

            push!(data, t1)

            # plot the control points?
            if ctrlPoints

                xP = [controlPoints[i][1] for i in eachindex(controlPoints)]
                yP = [controlPoints[i][2] for i in eachindex(controlPoints)]

                i = ["P$(i-1)" for i in eachindex(controlPoints)]

                t2 = PlotlyJS.scatter(;
                    x=xP, y=yP, mode="lines+markers+text", markersize=0.4, legend=:none, text=i, textposition="top center"
                )

                push!(data, t2)
            end

            PlotlyJS.plot(data)
        end

        @eval function plotSurface(S, controlPoints; ctrlPoints=true)

            data = PlotlyJS.GenericTrace[]

            x = [S[i, j][1] for i in eachindex(S[:, 1]), j in eachindex(S[1, :])]
            y = [S[i, j][2] for i in eachindex(S[:, 1]), j in eachindex(S[1, :])]
            z = [S[i, j][3] for i in eachindex(S[:, 1]), j in eachindex(S[1, :])]

            t1 = PlotlyJS.surface(; z=z, x=x, y=y, intensitymode="cell", colorscale="Viridis", opacity=0.75)

            push!(data, t1)

            # plot the control points?
            if ctrlPoints
                for ind in eachindex(controlPoints[:, 1])
                    aux = controlPoints[:, ind]

                    xP = [aux[i][1] for i in eachindex(aux)]
                    yP = [aux[i][2] for i in eachindex(aux)]
                    zP = [aux[i][3] for i in eachindex(aux)]

                    t2 = PlotlyJS.scatter3d(; x=xP, y=yP, z=zP, mode="markers+lines", markersize=0.4, legend=:none, showlegend=false)

                    push!(data, t2)
                end

                for ind in eachindex(controlPoints[1, :])
                    aux = controlPoints[ind, :]

                    xP = [aux[i][1] for i in eachindex(aux)]
                    yP = [aux[i][2] for i in eachindex(aux)]
                    zP = [aux[i][3] for i in eachindex(aux)]

                    t2 = PlotlyJS.scatter3d(; x=xP, y=yP, z=zP, mode="markers+lines", markersize=0.4, legend=:none, showlegend=false)

                    push!(data, t2)
                end
            end

            PlotlyJS.plot(data)
        end


    end
end
