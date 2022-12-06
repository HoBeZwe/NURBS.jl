
"""
    readMultipatch(filename::String)

Read multipatch file organized as:

PATCH a
4 4             -> degree in u and v
5 5             -> number of control points in u and v
0.0 0.0 ...     -> knot vector in u
0.0 0.0 ...     -> knot vector in v
0.1 0.2 ...     |
0.1 0.2 ...      > xyz components of the control points (normalized with the weigths -> we remove this weighting when reading in the data)
0.1 0.2 ...     |
1.0 1.0 ...     -> weights

PATCH b 
...
"""
function readMultipatch(filename::String, T=Float64)

    Patches = []

    # read the file into a vector of strings (one string per line)
    stringVec = readlines(filename)

    # loop over all lines
    lineInd = 1
    while lineInd ≤ length(stringVec)

        if startswith(stringVec[lineInd], "PATCH") # find PATCH and count from there

            Patch = parseSinglePatch(stringVec, lineInd, T)
            push!(Patches, Patch)

            lineInd += 9
            continue
        end

        lineInd += 1 # continue till next patch
    end

    return Patches
end


"""
    parseSinglePatch(stringVec, lineInd)

Take the string vector and the line index where "PATCH" stands and extract the information of a single patch
"""
function parseSinglePatch(stringVec, lineInd::Int, T)

    # --- read degrees
    degrees = parseLine(stringVec[lineInd + 1], Int)

    # --- read knot vectors
    uKnotVec = parseLine(stringVec[lineInd + 3], T)
    vKnotVec = parseLine(stringVec[lineInd + 4], T)

    # --- read control points
    pointDims = parseLine(stringVec[lineInd + 2], Int)
    controlPoints = parseCtrlPoints(pointDims, stringVec[(lineInd + 5):(lineInd + 7)], T)

    # --- read weights
    weights = parseLine(stringVec[lineInd + 8], T)
    weights = reshape(weights, pointDims[1], pointDims[2])

    # --- fill structure
    return NURBSsurface(Bspline(degrees[1], uKnotVec), Bspline(degrees[2], vKnotVec), controlPoints ./ weights, weights)
end


"""
    parseLine(line::String, T=Float64)

Take the string 'line' of the form '0.1 0.2 0.6 ...' and parse it to a vector with eltype 'T'.
"""
function parseLine(line::String, T=Float64)

    vector = T[]

    vecEntries = split(line, " "; keepempty=false) # split at blank spaces
    for (ind, vE) in enumerate(vecEntries) # parse each entry to a number and push it to the vector
        push!(vector, parse(T, vE))
    end

    return vector
end


"""
    parseCtrlPoints(pointDims, stringVec, T=Float64)

Take the 3 strings in 'stringVec' and parse it into a pointDims[1] x pointDims[2] matrix where each entry is a SVector for a controlpoint.
"""
function parseCtrlPoints(pointDims, stringVec, T=Float64)

    Nu = pointDims[1]
    Nv = pointDims[2]

    # read 3 lines containing the points
    auxX = parseLine(stringVec[1], T)
    auxY = parseLine(stringVec[2], T)
    auxZ = parseLine(stringVec[3], T)

    # sort into matrix of size (Nu, Nv)
    return [SVector(auxX[i + (j - 1) * Nu], auxY[i + (j - 1) * Nu], auxZ[i + (j - 1) * Nu]) for i in 1:Nu, j in 1:Nv]
end
