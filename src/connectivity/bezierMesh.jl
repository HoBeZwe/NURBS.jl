

function cellCart2Lin(cellU, cellV, pInd, cU, cV)

    return cellV + (cellU - 1) * cV + (pInd - 1) * cU * cV
end

function cellCart2LinTruncated(cellU, cellV, pInd, cU, cV)

    (cellU < 1 || cellU > cU || cellV < 1 || cellV > cV) && return -1 # requested cell is out of boundaries

    return cellCart2Lin(cellU, cellV, pInd, cU, cV)
end

function cellLin2Cart(cell, cU, cV)

    pInd, rem = divrem(cell - 1, cU * cV)
    cellU, cellV = divrem(rem, cV)

    return pInd + 1, cellU + 1, cellV + 1
end


# for a single patch cU ≠ cV is ok, for multipatch cU = cV is assumed
function bezierAdjacency(interfaces, commonVtxs, cU, cV, nP)

    patchAdj = BezierCellAdjacency{Int}[]

    # adjacency information for the Bezier cells internal to the patches
    for pInd in 1:nP
        for cellU in 1:cU
            for cellV in 1:cV

                v = vertAdjacencyPatch(cellU, cellV, pInd, cU, cV)
                e = edgeAdjacencyPatch(cellU, cellV, pInd, cU, cV)

                push!(patchAdj, BezierCellAdjacency(v, e))
            end
        end
    end

    # add adjacency information across interfaces
    for (ifcId, ifc) in enumerate(interfaces)

        p1, p2 = NURBS.patchIDs(ifc)
        l1, l2 = NURBS.localEdges(ifc)

        adjacencyInterface!(patchAdj, p1, p2, l1, l2, ifc.reverse, cU, cV)
        adjacencyInterface!(patchAdj, p2, p1, l2, l1, ifc.reverse, cU, cV)
    end

    # add adjacency information across touching corners of patches
    for (ifcId, ifc) in enumerate(commonVtxs)

        p1, p2 = NURBS.patchID(ifc[1]), NURBS.patchID(ifc[2])
        l1, l2 = NURBS.localEdge(ifc[1]), NURBS.localEdge(ifc[2])

        adjacencyCorner!(patchAdj, p1, p2, l1, l2, cU, cV)
        adjacencyCorner!(patchAdj, p2, p1, l2, l1, cU, cV)
    end

    return patchAdj
end


function adjacencyInterface!(patchAdj, p1Id, p2Id, l1, l2, rev, cU, cV)

    if l1 == 1
        for cellU in 1:cU
            cellInd = cellCart2Lin(cellU, 1, p1Id, cU, cV) # the cell to enter the information
            interfaceBezierCells!(patchAdj[cellInd], p2Id, cellU, l1, l2, rev, cU, cV, 1, 4)
        end

    elseif l1 == 2
        for cellV in 1:cV
            cellInd = cellCart2Lin(1, cellV, p1Id, cU, cV) # the cell to enter the information
            interfaceBezierCells!(patchAdj[cellInd], p2Id, cellV, l1, l2, rev, cU, cV, 1, 2)
        end

    elseif l1 == 3
        for cellU in 1:cU
            cellInd = cellCart2Lin(cellU, cV, p1Id, cU, cV) # the cell to enter the information
            interfaceBezierCells!(patchAdj[cellInd], p2Id, cellU, l1, l2, rev, cU, cV, 2, 3)
        end

    elseif l1 == 4
        for cellV in 1:cV
            cellInd = cellCart2Lin(cU, cellV, p1Id, cU, cV) # the cell to enter the information
            interfaceBezierCells!(patchAdj[cellInd], p2Id, cellV, l1, l2, rev, cU, cV, 4, 3)
        end
    end

    return nothing
end

function interfaceBezierCells!(patchAdj, p2Id, cell, l1, l2, rev, cU, cV, lE1, lE2)

    touchEdg = touchingEdge(p2Id, cell, l2, rev, cU, cV) # the touching cell at the edge
    tv1, tv2 = touchingVert(p2Id, cell, l2, rev, cU, cV) # the touching cells at the vertices

    push!(patchAdj.atLocalEdges[l1], touchEdg)

    tv1 > 0 && push!(patchAdj.atLocalVerts[lE1], tv1)
    tv2 > 0 && push!(patchAdj.atLocalVerts[lE2], tv2)

    return nothing
end


# determine the ajacent cell touching at an edge from the other patch of an interface
function touchingEdge(pID, cell, localEdge, reverse, cU, cV)

    i = reverseOrder(Val(reverse), cell, cV) # cV = cU assumed!

    localEdge == 1 && return cellCart2Lin(i, 1, pID, cU, cV)
    localEdge == 2 && return cellCart2Lin(1, i, pID, cU, cV)
    localEdge == 3 && return cellCart2Lin(i, cV, pID, cU, cV)
    localEdge == 4 && return cellCart2Lin(cU, i, pID, cU, cV)

    return -1 # localEdge ∉ {1,2,3,4}
end

reverseOrder(::Val{true}, cell, cV) = cV - cell + 1
reverseOrder(::Val{false}, cell, cV) = cell

# determine the ajacent cells touching at a vertex from the other patch of an interface
function touchingVert(pID, cell, localEdge, reverse, cU, cV)

    i = reverseOrder(Val(reverse), cell, cV) # cV = cU assumed!
    ip, is = offsetByOne(Val(reverse), i)

    localEdge == 1 && return cellCart2LinTruncated(is, 1, pID, cU, cV), cellCart2LinTruncated(ip, 1, pID, cU, cV)
    localEdge == 2 && return cellCart2LinTruncated(1, is, pID, cU, cV), cellCart2LinTruncated(1, ip, pID, cU, cV)
    localEdge == 3 && return cellCart2LinTruncated(is, cV, pID, cU, cV), cellCart2LinTruncated(ip, cV, pID, cU, cV)
    localEdge == 4 && return cellCart2LinTruncated(cU, is, pID, cU, cV), cellCart2LinTruncated(cU, ip, pID, cU, cV)

    return -1, -1 # localEdge ∉ {1,2,3,4}
end

offsetByOne(::Val{true}, i) = i - 1, i + 1
offsetByOne(::Val{false}, i) = i + 1, i - 1

function adjacencyCorner!(patchAdj, p1Id, p2Id, l1, l2, cU, cV)

    cellInd1 = cornerBezierCell(l1, p1Id, cU, cV) # the cell to enter the information
    cellInd2 = cornerBezierCell(l2, p2Id, cU, cV) # the touching cell

    push!(patchAdj[cellInd1].atLocalVerts[l1], cellInd2)

    return nothing
end

function cornerBezierCell(localEdge, pId, cU, cV)

    localEdge == 1 && return cellCart2Lin(1, 1, pId, cU, cV)
    localEdge == 2 && return cellCart2Lin(1, cV, pId, cU, cV)
    localEdge == 3 && return cellCart2Lin(cU, cV, pId, cU, cV)
    localEdge == 4 && return cellCart2Lin(cU, 1, pId, cU, cV)

    return -1 # localEdge ∉ {1,2,3,4}
end


# Determine for a cell the adjacent cells of the same patch at each local vertex 
function vertAdjacencyPatch(cellU, cellV, pInd, cU, cV)

    if cellU == 1 || cellV == 1
        v1 = Set{Int}()
    else
        v1 = Set(cellCart2Lin(cellU - 1, cellV - 1, pInd, cU, cV))
    end

    if cellU == 1 || cellV == cV
        v2 = Set{Int}()
    else
        v2 = Set(cellCart2Lin(cellU - 1, cellV + 1, pInd, cU, cV))
    end

    if cellU == cU || cellV == cV
        v3 = Set{Int}()
    else
        v3 = Set(cellCart2Lin(cellU + 1, cellV + 1, pInd, cU, cV))
    end

    if cellU == cU || cellV == 1
        v4 = Set{Int}()
    else
        v4 = Set(cellCart2Lin(cellU + 1, cellV - 1, pInd, cU, cV))
    end

    return [v1, v2, v3, v4]
end

# Determine for a cell the adjacent cells of the same patch at each local edge 
function edgeAdjacencyPatch(cellU, cellV, pInd, cU, cV)

    if cellV == 1
        e1 = Set{Int}()
    else
        e1 = Set(cellCart2Lin(cellU, cellV - 1, pInd, cU, cV))
    end

    if cellU == 1
        e2 = Set{Int}()
    else
        e2 = Set(cellCart2Lin(cellU - 1, cellV, pInd, cU, cV))
    end

    if cellV == cV
        e3 = Set{Int}()
    else
        e3 = Set(cellCart2Lin(cellU, cellV + 1, pInd, cU, cV))
    end

    if cellU == cU
        e4 = Set{Int}()
    else
        e4 = Set(cellCart2Lin(cellU + 1, cellV, pInd, cU, cV))
    end

    return [e1, e2, e3, e4]
end
