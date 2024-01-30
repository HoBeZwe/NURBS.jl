
# NOTE: have a closer look at opencascade for reading step files


"""
    readStep(filename::String, T=Float64)

Read a step file. 

So far only B_SPLINE_SURFACE_WITH_KNOTS and BOUNDED_SURFACE() are supported.
"""
function readStep(filename::String, T=Float64)

    # read the file into a vector of strings (one string per line)
    stringVec = readlines(filename)

    # extract Cartesian points
    points, offset = getCartesianPoints(stringVec, T)

    # initialize vector for the patches
    Patches = NURBSsurface{T}[] #Patches = BsplineSurface{T}[]

    # loop over all lines
    lineInd = 1
    while lineInd ≤ length(stringVec)

        # --- handle NURBS surfaces
        if startswith(stringVec[lineInd], "B_SPLINE_SURFACE") 

            lineCtlPts, lineInd = getCompleteStringTill(stringVec, lineInd, "B_SPLINE_SURFACE_WITH_KNOTS")
            lineKnVecs, lineInd = getCompleteStringTill(stringVec, lineInd, "GEOMETRIC_REPRESENTATION_ITEM")
            lineWeight, lineInd = getCompleteStringTill(stringVec, lineInd+1, "REPRESENTATION_ITEM")

            Patch = parse_NURBS_step(lineCtlPts, lineKnVecs, lineWeight, points, offset, T)
            push!(Patches, Patch)
        end

        # --- handle B-spline surfaces 
        if contains(stringVec[lineInd], "B_SPLINE_SURFACE_WITH_KNOTS") 

            line, lineInd = getCompleteString(stringVec, lineInd) # extract everything till ';' and return the reached line index
            Patch = parse_B_SPLINE_SURFACE_WITH_KNOTS(line, points, offset, T) # promoted to NURBS surface (weights are set to 1)
            push!(Patches, Patch)
        end

        lineInd += 1 # continue till next patch
    end

    return Patches
end


"""
parse_NURBS(lineCtlPts::String, lineKnVecs::String, lineWeight::String, points, offset::Int, T=Float64)

Example:

#163=(
BOUNDED_SURFACE()
B_SPLINE_SURFACE(4,4,((#349,#350,#351,#352,#353),(#354,#355,#356,#357,#358),
(#359,#360,#361,#362,#363),(#364,#365,#366,#367,#368),(#369,#370,#371,#372,
#373)),.UNSPECIFIED.,.F.,.F.,.F.)
B_SPLINE_SURFACE_WITH_KNOTS((5,5),(5,5),(0.,1.),(0.,1.),.UNSPECIFIED.)
GEOMETRIC_REPRESENTATION_ITEM()
RATIONAL_B_SPLINE_SURFACE(((5.07179676972449,4.52004210360334,4.35726558990816,
4.52004210360334,5.07179676972449),(4.52004210360334,3.86602540378444,3.64492370567392,
3.86602540378444,4.52004210360334),(4.35726558990816,3.64492370567392,3.40455735015306,
3.64492370567392,4.35726558990816),(4.52004210360334,3.86602540378444,3.64492370567392,
3.86602540378444,4.52004210360334),(5.07179676972449,4.52004210360334,4.35726558990816,
4.52004210360334,5.07179676972449)))
REPRESENTATION_ITEM('')
SURFACE()
);
"""
function parse_NURBS_step(lineCtlPts::String, lineKnVecs::String, lineWeights::String, points, offset::Int, T=Float64)

    firstBracketInd = findfirst('(', lineCtlPts) # find first comma

    # --- extract degrees 
    uDegree = parse(Int, lineCtlPts[firstBracketInd + 1]) # assume degree is only one digit
    vDegree = parse(Int, lineCtlPts[firstBracketInd + 3]) # assume degree is only one digit

    # --- extract control points
    controlPoints, pStopInd = extractCtrlPoints(lineCtlPts, firstBracketInd, points, offset, T)

    # --- extract knot vectors
    uKnotVec, vKnotVec = extractKnotVecs(lineKnVecs, 29, T) 

    # --- extract weigths
    pStartInd = findnext("(((", lineWeights, 25)[end]
    pStopInd  = findnext("))", lineWeights, 25)[1]

    # 5.071,4.520,4.357,4.520),(4.520,3.866,3.644,3.866),(4.520,3.866,3.644,3.866
    weightsString = lineWeights[(pStartInd + 1):(pStopInd - 1)]

    weights = ones(T, size(controlPoints))
    m = 1
    for list in eachsplit(weightsString, "),(")
        n = 1
        for entry in eachsplit(list, ',')

            weights[m,n] = parse(T, entry)
            n += 1
        end
        m += 1
    end

    return NURBSsurface(Bspline(uDegree, uKnotVec), Bspline(vDegree, vKnotVec), controlPoints, weights)
end


"""
    parse_B_SPLINE_SURFACE_WITH_KNOTS(line::String, points, offset::Int, T=Float64)

Parse a B_SPLINE_SURFACE_WITH_KNOTS

Example:

#13764=B_SPLINE_SURFACE_WITH_KNOTS('',3,3,((#16247,#16248,#16249,#16250,
#16251),(#16252,#16253,#16254,#16255,#16256),(#16257,#16258,#16259,#16260,
#16261),(#16262,#16263,#16264,#16265,#16266),(#16267,#16268,#16269,#16270,
#16271)),.UNSPECIFIED.,.F.,.F.,.F.,(4,1,4),(4,1,4),(0.,1.,2.),(0.,1.,2.),
 .UNSPECIFIED.);
"""
function parse_B_SPLINE_SURFACE_WITH_KNOTS(line::String, points, offset::Int, T=Float64)

    firstCommaInd = findfirst(',', line) # find first comma

    # --- extract degrees
    uDegree = parse(Int, line[firstCommaInd + 1])
    vDegree = parse(Int, line[firstCommaInd + 3])

    # --- extract control points
    controlPoints, pStopInd = extractCtrlPoints(line, firstCommaInd, points, offset, T)

    # --- extract knot vectors
    nextInd = pStopInd # --- skip next 4 entries
    for i in 1:5
        nextInd = findnext(',', line, nextInd) + 1
    end

    uKnotVec, vKnotVec = extractKnotVecs(line, nextInd, T)    

    #return BsplineSurface(Bspline(uDegree, uKnotVec), Bspline(vDegree, vKnotVec), controlPoints) 
    return NURBSsurface(Bspline(uDegree, uKnotVec), Bspline(vDegree, vKnotVec), controlPoints, ones(size(controlPoints)))
end


"""
    getCartesianPoints(stringVec)

Extract all Cartesian points and the offset between the first point and the #counting.
(Assumption: all points appear in a consecutive list.)
"""
function getCartesianPoints(stringVec::Vector{String}, T=Float64)

    # loop over all lines
    lineInd = 1
    offset = -1

    points = []
    while lineInd ≤ length(stringVec)

        line = stringVec[lineInd]

        if contains(line, "CARTESIAN_POINT") # find CARTESIAN_POINT

            # determine difference between line number and #counter
            if offset < 0
                endInd = findfirst('=', line) - 1
                offset = parse(Int, line[2:endInd]) - 1
            end

            # extract point
            nextInd = findfirst(',', line) + 1
            point, aux = extractBracketList(line, nextInd, T)
            push!(points, point)
        end

        lineInd += 1 # continue till next patch
    end

    return points, offset
end


"""
    getCompleteString(stringVec::String, nextLineInd::Int)

Retreive all lines of the string vector 'stringVec' starting from 'nextLineInd' till a semicolon is encountered.
"""
function getCompleteString(stringVec::Vector{String}, nextLineInd::Int)

    completeString = ""

    while !endswith(stringVec[nextLineInd], ';')

        completeString = completeString * stringVec[nextLineInd]
        nextLineInd += 1
    end

    return completeString * stringVec[nextLineInd], nextLineInd
end


"""
    getCompleteStringTill(stringVec, nextLineInd::Int, stopWhenEncountered::String)

Retreive all lines of the string vector 'stringVec' starting from 'nextLineInd' till the next line starts with 'stopWhenEncountered'.
"""
function getCompleteStringTill(stringVec, nextLineInd::Int, stopWhenEncountered::String)

    completeString = ""

    while !startswith(stringVec[nextLineInd], stopWhenEncountered)

        completeString = completeString * stringVec[nextLineInd]
        nextLineInd += 1
    end

    return completeString, nextLineInd
end


"""
    extractBracketList(line::String, nextInd::Int, T=Float64)

Given a string with a list of numbers between two brackets (e.g., (2,3,1,5,...,6)), extract the numbers.

As second output the index after ')' is returned.

'nextInd' is the index of the '(' in 'line'.
"""
function extractBracketList(line::String, nextInd::Int, T=Float64)

    indStart = findnext('(', line, nextInd) + 1
    indStop  = findnext(')', line, nextInd) - 1

    list = T[]
    for entry in eachsplit(line[indStart:indStop], ',')
        push!(list, parse(T, entry))
    end

    return list, indStop + 2
end


"""
    extractCtrlPoints(line::String, firstCommaInd::Int, points, offset::Int, T=Float64)

Extract a list of points between '((' and '))' in the line starting from 'firstCommaInd' of the format 

        ((#349,#350,#351,#352,#353),(#354,#355,#356,#357,#358),(#359,#360,#361,#362,#363),(#364,#365,#366,#367,#368),(#369,#370,#371,#372,#373))
"""
function extractCtrlPoints(line::String, firstCommaInd::Int, points, offset::Int, T=Float64)

    pStartInd = findnext("((", line, firstCommaInd)[end]
    pStopInd  = findnext("))", line, firstCommaInd)[1]

    # #766,#767,#768,#769,#770),(#771,#772,#773,#774,#775),(#776,#777,#778,#779,#780),(#781,#782,#783,#784,#785),(#786,#787,#788,#789,#790
    ctrlPointsIndsString = line[(pStartInd + 1):(pStopInd - 1)]


    ctrlPointsInds = []
    for list in eachsplit(ctrlPointsIndsString, "),(")

        listRem = replace(list, "#" => "") # remove # symbols

        ctrPointsInd = Int[]
        for entry in eachsplit(listRem, ",")
            push!(ctrPointsInd, parse(Int, entry))
        end

        push!(ctrlPointsInds, ctrPointsInd)
    end

    sizeV = length(ctrlPointsInds)
    sizeU = length(ctrlPointsInds[1])

    return [SVector{3,T}(points[ctrlPointsInds[v][u] - offset]) for v in 1:sizeV, u in 1:sizeU], pStopInd
end


"""
    extractKnotVecs(line::String, startInd::Int, T=Float64)

Extract the knot vectors from the 'line' of the format 'SOMESTRING(5,5),(5,5),(0.,1.),(0.,1.)SOMESTRING' where 'startInd' points to the 'G'.
"""
function extractKnotVecs(line::String, startInd::Int, T=Float64)

    uMultipli, nextInd = extractBracketList(line, startInd, Int)
    vMultipli, nextInd = extractBracketList(line, nextInd, Int)

    uKvec, nextInd = extractBracketList(line, nextInd, T)
    vKvec, nextInd = extractBracketList(line, nextInd, T)

    uKnotVec = T[]
    for i in eachindex(uKvec)
        append!(uKnotVec, repeat([uKvec[i]], uMultipli[i])) # repeat according to multiplicity
    end

    vKnotVec = T[]
    for i in eachindex(vKvec)
        append!(vKnotVec, repeat([vKvec[i]], vMultipli[i])) # repeat according to multiplicity
    end

    return uKnotVec, vKnotVec
end