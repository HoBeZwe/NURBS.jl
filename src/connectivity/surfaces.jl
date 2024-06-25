
"""

Local edge and vertex numbering

         2   2
  1 o--------o-----> v-direction
    |        |
   1|        |3
    |        |
  4 o--------o 3
    |    4
    v
    u-direction

An interface of a patch.
"""
struct PatchInterface{I}
    patchID::I
    localEdge::I
end


"""

An interface is defined by two patches: their ID and wether the vectors normal to the interface point in the same direction.
"""
struct Interface{I}
    patch1::PatchInterface{I}
    patch2::PatchInterface{I}
    orientation::I              # do basis functions normal to the edge point in same direction (+1 for yes, -1 for no)
    reverse::Bool               # are the basis functions on the two patches sorted in the same order
    globalVIds::SVector{2,I}
end

patchID(pifc::PatchInterface) = pifc.patchID
localEdge(pifc::PatchInterface) = pifc.localEdge

patchIDs(ifc::Interface) = patchID(ifc.patch1), patchID(ifc.patch2)
localEdges(ifc::Interface) = localEdge(ifc.patch1), localEdge(ifc.patch2)

"""
    InterfacePatchwise

Save connectivity information for a single patch: 

Each entry is a vector, e.g., localEdgesVec[i] contains the local edges being part of an interface of patch i.
"""
struct InterfacePatchwise
    localEdgesVec   # the local edges being part of an interface
    orientationVec  # the orientations of the interface functions with respect to the functions on the patch (±1)
    interfaceIdVec  # the Id's of the interfaces
    reversionVec    # do the basis functions have to be sorted in reverse order (bool)
    localVtxsVec    # local vertices of the patch where the patch touches another patch at the corner only
    commVtxIdVec    # the Id's of the vertices where the patch touches another patch at the corner only
end

struct BezierCellAdjacency{I}
    atLocalVerts::Vector{Set{I}}
    atLocalEdges::Vector{Set{I}}
end


struct InterfaceData{I}
    ifc::Vector{Interface{I}}
    patchIfs::InterfacePatchwise
end

struct Connectivity{I}
    interfaces::Vector{Interface{I}}
    patchInterfaces::InterfacePatchwise
    bezierAdjacency::Vector{BezierCellAdjacency{I}}
end

function Connectivity(patches::Vector{<: Surface}, cU, cV; tol=1e-3) 

    interfaces, commonVtxs = identifyInterfaces(patches, tol=tol)
    ifPatchwise = getPatchInterfaces(patches, interfaces, commonVtxs)
    bezierAdj = bezierAdjacency(interfaces, commonVtxs, cU, cV, length(patches))

    return Connectivity(interfaces, ifPatchwise, bezierAdj)
end


"""

Determine the interfaces between patches including the property if the basis functions normal to the edge point in the same direction.

Note: not an efficient implementation. Preferably write the information in the geometry file.
"""
function identifyInterfaces(patches::Vector{<:Surface{F}}; tol=1e-3) where {F}

    T = Int
    Intrfcs = Interface{T}[] # array to store the interfaces
    soFar = Set{Vector{T}}() # auxiliary array to keep track of which interfaces where already found

    x = T(0)
    commonVtxs = typeof((PatchInterface(x, x), PatchInterface(x, x)))[]
    soFarVtxs = Set{Vector{T}}()

    uniqueVerts = uniqueVertices(patches)

    for (i, P) in enumerate(patches) # loop over all patches

        Pctrl = P.controlPoints

        P11 = Pctrl[1, 1]
        P21 = Pctrl[end, 1]
        P12 = Pctrl[1, end]
        P22 = Pctrl[end, end]

        for (j, Pother) in enumerate(patches) # loop over all patches

            i == j && continue # do not compare patch with itself

            PotherCtrl = Pother.controlPoints

            P11_ = PotherCtrl[1, 1] # only compare corner points (maybe generalize in future)
            P21_ = PotherCtrl[end, 1]
            P12_ = PotherCtrl[1, end]
            P22_ = PotherCtrl[end, end]

            # do patches have two points in common? (always either 0, 1, or 2 everything else would be a non-manifold)
            pInds = T[]
            pIndsOther = T[]

            for (inda, Pa) in enumerate([P11, P12, P22, P21])
                for (indb, Pb) in enumerate([P11_, P12_, P22_, P21_])

                    if norm(Pa - Pb) < tol # points agree

                        push!(pInds, inda)
                        push!(pIndsOther, indb)
                    end
                end
            end

            nCommonPts = length(pInds)

            if nCommonPts == 0 # patches do not touch

                continue

            elseif nCommonPts == 1 # patches have common vertex

                # check whether the common vertex was already found
                alreadyOccurred = false
                for (ind, ss) in enumerate(soFarVtxs)
                    alreadyOccurred |= issetequal([i, j], ss) # compare with the already found interfaces
                end

                if !alreadyOccurred
                    push!(soFarVtxs, [i, j])
                    push!(commonVtxs, (PatchInterface(i, pInds[1]), PatchInterface(j, pIndsOther[1])))
                end

            elseif nCommonPts == 2 # patches have common edge

                issetequal(pInds, [1, 2]) &&
                    (edge = 2; glVerts = [getGlobalIndex(P11, uniqueVerts; tol=tol), getGlobalIndex(P12, uniqueVerts; tol=tol)])
                issetequal(pInds, [2, 3]) &&
                    (edge = 3; glVerts = [getGlobalIndex(P12, uniqueVerts; tol=tol), getGlobalIndex(P22, uniqueVerts; tol=tol)])
                issetequal(pInds, [3, 4]) &&
                    (edge = 4; glVerts = [getGlobalIndex(P22, uniqueVerts; tol=tol), getGlobalIndex(P21, uniqueVerts; tol=tol)])
                issetequal(pInds, [4, 1]) &&
                    (edge = 1; glVerts = [getGlobalIndex(P21, uniqueVerts; tol=tol), getGlobalIndex(P11, uniqueVerts; tol=tol)])

                issetequal(pIndsOther, [1, 2]) && (otherEdge = 2)
                issetequal(pIndsOther, [2, 3]) && (otherEdge = 3)
                issetequal(pIndsOther, [3, 4]) && (otherEdge = 4)
                issetequal(pIndsOther, [4, 1]) && (otherEdge = 1)

                getGlobalIndex(P22, uniqueVerts; tol=tol)

                # basis functions normal to edge point in same direction on both patches for combinations 11, 22, 33, 44, 12, 21, 34, 43 denoting the local edge numbers in the adjacent patches
                if (edge == otherEdge || issetequal([edge, otherEdge], [1, 2]) || issetequal([edge, otherEdge], [3, 4]))
                    orientation = -1
                else # 13, 31, 14, 41, 23, 32, 24, 42
                    orientation = +1
                end

                # basis functions normal to edge are sorted in oposite order on the patches for combinations 11, 22, 33, 44, 32, 23, 14, 41 denoting the local edge numbers in the adjacent patches
                if (edge == otherEdge || issetequal([edge, otherEdge], [2, 3]) || issetequal([edge, otherEdge], [1, 4]))
                    reverse = true
                else # 13, 31, 14, 41, 23, 32, 24, 42
                    reverse = false
                end

                # check whether the interface was already found
                alreadyOccurred = false
                for (ind, ss) in enumerate(soFar)
                    alreadyOccurred |= issetequal([i, j], ss) # compare with the already found interfaces
                end

                if !alreadyOccurred
                    push!(soFar, [i, j])
                    push!(
                        Intrfcs,
                        Interface(PatchInterface(i, edge), PatchInterface(j, otherEdge), orientation, reverse, SVector{2}(glVerts)),
                    )
                end

            else
                error("Interfaces intersect!\n")
            end
        end
    end

    return Intrfcs, commonVtxs
end


"""


Determine for each patch the local edges where there is an interface and whether the orientation has to be flipped.

    localEdgesVec[p] contains the local edges of the p-th patch that are attached to an interface
"""
function getPatchInterfaces(patches::Vector{<:Surface{F}}, interfaces, commonVtxs) where {F}

    localEdgesVec = Vector{Vector{Int}}()
    orientationVec = Vector{Vector{Int}}()
    interfaceIdVec = Vector{Vector{Int}}()
    reversionVec = Vector{Vector{Bool}}()

    localVtxsVec = Vector{Vector{Int}}()
    commVtxIdVec = Vector{Vector{Int}}()

    for patchID in eachindex(patches)

        localEdges = Int[]
        orientation = Int[]
        interfaceId = Int[]
        reversion = Bool[]
        localVtxs = Int[]
        commVtxId = Int[]

        for (ind, b) in enumerate(interfaces)

            if b.patch1.patchID == patchID
                push!(localEdges, b.patch1.localEdge)
                push!(orientation, 1)
                push!(interfaceId, ind)
                push!(reversion, b.reverse)

            elseif b.patch2.patchID == patchID
                push!(localEdges, b.patch2.localEdge)
                push!(orientation, b.orientation)
                push!(interfaceId, ind)
                push!(reversion, b.reverse)
            end
        end

        for (ind, v) in enumerate(commonVtxs)

            if v[1].patchID == patchID
                push!(localVtxs, v[1].localEdge)
                push!(commVtxId, ind)

            elseif v[2].patchID == patchID
                push!(localVtxs, v[2].localEdge)
                push!(commVtxId, ind)
            end
        end

        push!(localEdgesVec, localEdges)
        push!(orientationVec, orientation)
        push!(interfaceIdVec, interfaceId)
        push!(reversionVec, reversion)

        push!(localVtxsVec, localVtxs)
        push!(commVtxIdVec, commVtxId)
    end

    return InterfacePatchwise(localEdgesVec, orientationVec, interfaceIdVec, reversionVec, localVtxsVec, commVtxIdVec)
end


"""
    uniqueVertices(Patches::Vector{<:Surface{F}}; tol=1e-3) where {F}

Retreive the unique corner points of the patches in the provided list.
"""
function uniqueVertices(Patches::Vector{<:Surface{F}}; tol=1e-3) where {F}

    Pctrl = Patches[1].controlPoints

    P11 = Pctrl[1, 1]
    P21 = Pctrl[end, 1]
    P12 = Pctrl[1, end]
    P22 = Pctrl[end, end]


    pointList = [P11, P21, P12, P22]

    for (pInd, patch) in enumerate(Patches)

        Pctrl = patch.controlPoints

        P11 = Pctrl[1, 1]
        P21 = Pctrl[end, 1]
        P12 = Pctrl[1, end]
        P22 = Pctrl[end, end]

        !isin(P11, pointList; tol=tol) && push!(pointList, P11)
        !isin(P21, pointList; tol=tol) && push!(pointList, P21)
        !isin(P12, pointList; tol=tol) && push!(pointList, P12)
        !isin(P22, pointList; tol=tol) && push!(pointList, P22)
    end

    return pointList
end


"""
    getGlobalIndex(P, Plist; tol=1e-3)

Get the global index of the point P with respect to the ordering in the list of points Plist.
"""
function getGlobalIndex(P, Plist; tol=1e-3)

    for (ind, Pentry) in enumerate(Plist)

        norm(Pentry - P) < tol && return ind
    end

    return error("point not found")
end


"""
    isin(P, Plist; tol=1e-3)

Check if 3D point P is in list of points Plist.
"""
function isin(P, Plist; tol=1e-3)

    for (ind, Pentry) in enumerate(Plist)

        if norm(P - Pentry) < tol # points agree
            return true
        end
    end

    return false
end
