
function coordis(S)
    x = [S[i, j][1] for i in eachindex(S[:, 1]), j in eachindex(S[1, :])]
    y = [S[i, j][2] for i in eachindex(S[:, 1]), j in eachindex(S[1, :])]
    z = [S[i, j][3] for i in eachindex(S[:, 1]), j in eachindex(S[1, :])]

    return x[:], y[:], z[:]
end


"""
    vtk(patches, resolution=0.01)

Prepare the patches to be displayed by Paraview: 

The output can be written to the VTK format using the WriteVTK.jl package.

using WriteVTK

vtk_grid("path/filename", x, y, z, cellV) do vtk
end
"""
function vtk(patches, resolution=0.01)

    uEvalpoints = collect(0:resolution:1.0)
    vEvalpoints = collect(0:resolution:1.0)

    Nv = length(vEvalpoints)
    Nu = length(uEvalpoints)
    Nc = length(patches)

    cellV = typeof(MeshCell(VTKCellTypes.VTK_QUAD, [2, 4, 3, 5]))[]
    for cell in 0:(Nc - 1)
        for j in 0:(Nu - 2), i in 1:(Nv - 1)
            ind = j * Nv + i + cell * Nu * Nv
            push!(cellV, MeshCell(VTKCellTypes.VTK_QUAD, [ind + Nv, ind + Nv + 1, ind + 1, ind]))
        end
    end

    x = Vector{Float64}()
    y = Vector{Float64}()
    z = Vector{Float64}()

    for (_, patch) in enumerate(patches)

        S = patch(uEvalpoints, vEvalpoints)'
        xs, ys, zs = coordis(S)

        append!(x, xs)
        append!(y, ys)
        append!(z, zs)
    end

    return cellV, x, y, z
end


"""
    saveVtk(filename::String, patches; resolution=0.01)

Save as an unstructured grid in a .vtu file to be displayed by ParaView.
"""
function saveVtk(filename::String, patches; resolution=0.01)

    cellV, x, y, z = vtk(patches, resolution)

    vtk_grid(filename, x, y, z, cellV) do vtk
    end

    return nothing
end
