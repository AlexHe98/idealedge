"""
Find the ideal edges after crushing a normal surface.
"""
from sys import argv
from regina import *
from loop import NotLoop, IdealLoop


#TODO
#_quadMapping = [
#        Perm4(0,1,2,3),
#        Perm4(0,2,1,3),
#        Perm4(0,3,1,2) ]
_quadSameSide = [
        Perm4(1,0,3,2),
        Perm4(2,3,0,1),
        Perm4(3,2,1,0) ]
_quadOpposite = [
        Perm4(2,3,0,1),
        Perm4(1,0,3,2),
        Perm4(1,0,3,2)]


def findIdealEdges( surf, targets=None ):
    """
    Returns details of the ideal edges after crushing surf.

    Specifically, this routine returns a list of pairs of the form (t, n),
    where t is the index after crushing of a tetrahedron that will be
    incident to some ideal edge E, and n is an edge number of this
    tetrahedron that corresponds to the edge E. Each ideal edge will be
    represented by exactly one such pair in the returned list.

    If the dictionary of surviving segments has been precomputed using the
    _survivingSegments() routine, then this can be supplied using the
    optional targets argument. Otherwise, this routine will compute the
    surviving segments for itself.
    """
    tri = surf.triangulation()
    #TODO Decide if we need targets.
    if targets is None:
        targets = _survivingSegments(surf)

    #NOTE BEGIN TRAVERSE VERTICAL
    # Traverse all vertical boundary components of the parallelity bundle.
    pairedVertices = [
            [ 1,0,3,2 ],
            [ 2,3,0,1 ],
            [ 3,2,1,0 ] ]
    alreadyTraversedVertBdryQuads = set()
    for tet in tri.tetrahedra():
        quadType = _quadType( surf, tet.index() )
        if quadType is None:
            continue

        # This tet contains quads, so it might touch...
        for face in range(3):
            if ( tet.index(), face ) in alreadyTraversedVertBdryQuads:
                continue

            # Found a vertical boundary component that we have yet to
            # traverse.
    #NOTE END TRAVERSE VERTICAL

    #NOTE BEGIN TRAVERSE SURVIVING
    #NOTE END TRAVERSE SURVIVING

    #TODO
    raise NotImplementedError()


def parallelityBundleWeights( surf, weights ):
    """
    Calculates the total weight of each component of the parallelity bundle.

    In detail, this routine returns a list consisting of pairs (s,w), where
    s is a list of type-2 segments that are all in the same component of the
    parallelity bundle, and w is the total weight of that component.

    WARNING:
    --> This routine does not check whether the given segment weights make
        sense for the given normal surface. Therefore, this routine may
        produce undefined behaviour when given invalid input arguments.
    """
    # Naively compute orbits using repeated depth-first search.
    # At least in theory, it would be better to do this using something like
    # the Agol-Hass-Thurston weighted orbit-counting algorithm.
    segments = _segments(surf)
    orbitWeights = []
    while segments:
        seg, segType = segments.popitem()
        if segType != 2:
            continue

        # Compute the total weight of the orbit of this type-2 segment.
        orbit = _computeOrbit( surf, seg )
        totalWeight = 0
        for seg in orbit:
            totalWeight += weights[seg]
        orbitWeights.append( ( orbit, weight ) )
    return orbitWeights


def _computeOrbit( surf, start ):
    tri = surf.triangulation()  #TODO Is this needed?

    # Depth-first search.
    stack = [start]
    visited = set()
    while stack:
        current = stack.pop()
        if current in visited:
            continue

        # We haven't visited the current segment yet, so we need to find all
        # segments that are adjacent to it along parallel cells or faces.
        #TODO Update implementation and documentation.
        ei, seg = current
        e = tri.edge(ei)
        wt = surf.edgeWeight(ei).safeLongValue()
        visited.add(current)    # Record as visited now, so we don't forget.
        for emb in e.embeddings():
            tet = emb.tetrahedron()
            teti = tet.index()
            en = emb.face()
            ver = emb.vertices()

            # To locate the relevant parallel cells and faces in tet, need to
            # get the normal coordinates incident to e.
            f = [ surf.triangles( teti, ver[i] ).safeLongValue()
                    for i in range(2) ]
            q = 0
            qType = None
            for qt in range(3):
                if qt in { en, 5-en }:
                    continue
                quads = surf.quads( teti, qt ).safeLongValue()
                if quads > 0:
                    q = quads
                    qType = qt
                    break

            # Does the current segment belong to a parallel cell or face in
            # this tet? If so, then we need to find all adjacent segments.
            if seg < f[0]:
                # The current segment belongs to a parallel triangular cell
                # at vertex ver[0].
                for otherEnd in range(4):
                    if otherEnd in { ver[0], ver[1] }:
                        continue

                    # The current segment is adjacent to a segment of the
                    # edge with endpoints ver[0] and otherEnd.
                    eiOther = tet.edge( ver[0], otherEnd ).index()
                    enOther = Edge3.edgeNumber[ver[0]][otherEnd]
                    verOther = tet.edgeMapping(enOther)
                    if verOther[0] == ver[0]:
                        adjacent = ( eiOther, seg )
                    else:
                        wtOther = surf.edgeWeight(eiOther).safeLongValue()
                        adjacent = ( eiOther, wtOther - seg )

                    # If the adjacent segment is one of the targets, then we
                    # are done; otherwise, we add it to the stack.
                    output = targets.get( adjacent, None )
                    if output is not None:
                        found.add(output)
                        #return output
                    else:
                        stack.append(adjacent)
            elif seg > f[0] + q:
                # The current segment belongs to a parallel triangular cell
                # at vertex ver[1].
                for otherEnd in range(4):
                    if otherEnd in { ver[0], ver[1] }:
                        continue

                    # The current segment is adjacent to a segment of the
                    # edge with endpoints ver[1] and otherEnd.
                    eiOther = tet.edge( ver[1], otherEnd ).index()
                    enOther = Edge3.edgeNumber[ver[1]][otherEnd]
                    verOther = tet.edgeMapping(enOther)
                    if verOther[0] == ver[1]:
                        adjacent = ( eiOther, wt - seg )
                    else:
                        wtOther = surf.edgeWeight(eiOther).safeLongValue()
                        adjacent = ( eiOther, wtOther - wt + seg )

                    # If the adjacent segment is one of the targets, then we
                    # are done; otherwise, we add it to the stack.
                    output = targets.get( adjacent, None )
                    if output is not None:
                        found.add(output)
                        #return output
                    else:
                        stack.append(adjacent)
            elif q > 0:
                # The quadrilaterals divide tet into two "sides". The edge
                # opposite this segment has endpoints lying on different
                # sides, so we label these endpoints accordingly.
                side = [ { 0, qType + 1 } ]
                side.append( {0,1,2,3} - side[0] )
                if ver[0] not in side[0]:
                    side[0], side[1] = side[1], side[0]
                side[0].remove( ver[0] )
                side[1].remove( ver[1] )
                opp = [ side[0].pop(), side[1].pop() ]

                # Find all edges containing segments that are adjacent to the
                # current segment.
                qDepth = seg - f[0]     # 0 <= qDepth <= q
                if qDepth == 0:
                    adjEndpoints = [ [ ver[0], opp[1] ] ]
                elif qDepth == q:
                    adjEndpoints = [ [ opp[0], ver[1] ] ]
                else:
                    adjEndpoints = [
                            [ ver[0], opp[1] ],
                            [ opp[0], ver[1] ],
                            [ opp[0], opp[1] ] ]
                for start, end in adjEndpoints:
                    # The current segment is adjacent to a segment of the
                    # edge going from start to end.
                    eiAdj = tet.edge( start, end ).index()
                    enAdj = Edge3.edgeNumber[start][end]
                    verAdj = tet.edgeMapping(enAdj)
                    triangles = surf.triangles(
                            teti, verAdj[0] ).safeLongValue()
                    if verAdj[0] == start:
                        adjacent = ( eiAdj, triangles + qDepth )
                    else:
                        adjacent = ( eiAdj, triangles + q - qDepth )

                    # If the adjacent segment is one of the targets, then we
                    # are done; otherwise, we add it to the stack.
                    output = targets.get( adjacent, None )
                    if output is not None:
                        found.add(output)
                        #return output
                    else:
                        stack.append(adjacent)

        # End of loop. Move on to the next embedding of edge e.

    # If the search terminates without finding the ideal edge, then the ideal
    # edge must belong to a component that gets destroyed.
    #TODO
    raise NotImplementedError()
    #return found
    #return None


def _segments(surf):
    """
    Returns all segments induced by the given normal surface.

    In detail, the returned object is a dictionary whose items are of the
    form
        (e,s): t
    where:
    --> e is the index of the edge containing the segment;
    --> s is position of the segment along edge e; and
    --> t is the type of the segment.
    """
    tri = surf.triangulation()
    segments = dict()
    for ei in range( tri.countEdges() ):
        edgeWt = surf.edgeWeight(ei).safeLongValue()
        if edgeWt == 0:
            segments[ (ei,0) ] = 0
        else:
            for seg in range( edgeWt + 1 ):
                if seg in {0,edgeWt}:
                    segments[ (ei,seg) ] = 1
                else:
                    segments[ (ei,seg) ] = 2
    return segments


def _eulerWeights(surf):
    """
    Returns a dictionary that assigns weights to the edge segments induced by
    the given normal surface, such that the total weight in each component C
    of the parallelity bundle gives twice the Euler characteristic of the
    base of C.
    """
    tri = surf.triangulation()

    # Vertices of the base correspond to type-2 segments.
    weights = dict()
    segments = _segments(surf)
    for seg in segments:
        segType = segments[seg]
        if segType == 2:
            weights[seg] = 2
        else:
            #TODO TEST
            pass
            #weights[seg] = 0

    # Edges of the base correspond to parallel faces of the induced cell
    # decomposition, and faces of the base correspond to parallel cells.
    #
    # Let:  E = number of edges of the base
    #       F = number of faces of the base
    #       q = number of parallel quad cells
    #       t = number of parallel triangle cells
    #       b = number of parallel faces in the vertical boundary
    #       i = number of isolated parallel faces 
    # We adjust weights using the following observation:
    #   -2E + 2F = -(4q + 3t + b + 2i) + 2(q + t) = -2q - t - b - 2i.
    for tet in tri.tetrahedra():
        teti = tet.index()

        # For each parallel triangle cell, adjust the weight on one incident
        # segment.
        for triType in range(4):
            triCount = surf.triangles( teti, triType ).safeLongValue()
            if triCount < 2:
                continue

            # We have at least one parallel triangle cell.
            if triType in {0,3}:
                incidentEdgeNum = 2
            else:
                incidentEdgeNum = 3
            ei = tet.edge(incidentEdgeNum).index()
            if tet.edgeMapping(incidentEdgeNum)[0] == triType:
                for segment in range( 1, triCount ):
                    weights[ (ei,segment) ] -= 1
                    #TODO TEST
                    #print(teti,"tri -1")
            else:
                edgeWt = surf.edgeWeight(ei).safeLongValue()
                for segment in range( edgeWt - 1, edgeWt - triCount, -1 ):
                    weights[ (ei,segment) ] -= 1
                    #TODO TEST
                    #print(teti,"tri -1")

        # We also need to adjust weights using parallel quad cells, and
        # parallel faces that are either boundary or isolated. Such cases
        # only arise in tetrahedra that contain at least one quad.
        quads = _quads( surf, teti )
        if quads is None:
            continue
        quadType, quadCount = quads

        # For each parallel quad cell, adjust the weight on one incident
        # segment.
        if quadCount >= 2:
            # We have at least one parallel quad cell.
            incidentEdgeNum = quadType+2    # Can check this works by hand.
            ei = tet.edge(incidentEdgeNum).index()
            v = tet.edgeMapping(incidentEdgeNum)[0]
            quadStart = surf.triangles( teti, v ).safeLongValue()
            for segment in range( quadStart + 1, quadStart + quadCount ):
                weights[ (ei,segment) ] -= 2
                #TODO TEST
                #print(teti,"quad -2")

        # For each parallel face that is either isolated or on the boundary
        # of the parallelity bundle, adjust the weight on one incident
        # segment. In the isolated case, we only need to subtract 1 on this
        # side, because we will eventually subtract another 1 on the other
        # side as well.
        sameSide = _quadSameSide[quadType]
        opposite = _quadOpposite[quadType]
        for face in range(4):
            triCount = surf.triangles( teti, sameSide[face] ).safeLongValue()
            if triCount == 0:
                continue

            # We have a parallel face that is either isolated or boundary.
            incidentEdgeNum = Edge3.faceNumber( Perm4(
                sameSide[face], opposite[face],
                face, opposite[sameSide[face]] ) )
            ei = tet.edge(incidentEdgeNum).index()
            if tet.edgeMapping(incidentEdgeNum)[0] == sameSide[face]:
                segment = triCount
            else:
                edgeWt = surf.edgeWeight(ei).safeLongValue()
                segment = edgeWt - triCount
            weights[ (ei,segment) ] -= 1
            #TODO TEST
            #print(teti,"face -1")

    # All done!
    return weights


def _quadType( surf, tetIndex ):
    """
    Returns the quad type in which the given normal surface intersects the
    tetrahedron with the given index, or None if there is no such quad.
    """
    quads = _quads( surf, tetIndex )
    if quads is None:
        return None
    else:
        return quads[0]


def _quads( surf, tetIndex ):
    """
    Returns the quad type and the number of quads in which the given normal
    surface intersects the tetrahedron with the given index, or None if there
    is no such quad.
    """
    for quadType in range(3):
        quadCount = surf.quads( tetIndex, quadType ).safeLongValue()
        if quadCount > 0:
            return ( quadType, quadCount )
    return None


def findIdealEdges_old( surf, start, targets=None ):
    """
    Returns details of the ideal edge that corresponds to the given start
    segment after crushing surf.

    Specifically, if the ideal edge belongs to a component that gets
    destroyed after crushing, then this routine returns None. Otherwise, this
    routine returns a pair (t, e), where t is the index after crushing of a
    tetrahedron that will be incident to the ideal edge, and e is an edge
    number of this tetrahedron that corresponds to the ideal edge.

    If the dictionary of surviving segments has been precomputed using the
    _survivingSegments() routine, then this can be supplied using the
    optional targets argument. Otherwise, this routine will compute the
    surviving segments for itself.
    """
    #TODO Fix documentation/comments.
    tri = surf.triangulation()
    if targets is None:
        targets = _survivingSegments(surf)

    # If the start segment is one of the targets, then we are already done.
    found = set()
    output = targets.get( start, None )
    if output is not None:
        found.add(output)
        #return output

    # Otherwise, we find the ideal edge using depth-first search.
    stack = [start]
    visited = set()
    while stack:
        current = stack.pop()
        if current in visited:
            continue

        # We haven't visited the current segment yet, so we need to find all
        # segments that are adjacent to it along parallel cells or faces.
        ei, seg = current
        e = tri.edge(ei)
        wt = surf.edgeWeight(ei).safeLongValue()
        visited.add(current)    # Record as visited now, so we don't forget.
        for emb in e.embeddings():
            tet = emb.tetrahedron()
            teti = tet.index()
            en = emb.face()
            ver = emb.vertices()

            # To locate the relevant parallel cells and faces in tet, need to
            # get the normal coordinates incident to e.
            f = [ surf.triangles( teti, ver[i] ).safeLongValue()
                    for i in range(2) ]
            q = 0
            qType = None
            for qt in range(3):
                if qt in { en, 5-en }:
                    continue
                quads = surf.quads( teti, qt ).safeLongValue()
                if quads > 0:
                    q = quads
                    qType = qt
                    break

            # Does the current segment belong to a parallel cell or face in
            # this tet? If so, then we need to find all adjacent segments.
            if seg < f[0]:
                # The current segment belongs to a parallel triangular cell
                # at vertex ver[0].
                for otherEnd in range(4):
                    if otherEnd in { ver[0], ver[1] }:
                        continue

                    # The current segment is adjacent to a segment of the
                    # edge with endpoints ver[0] and otherEnd.
                    eiOther = tet.edge( ver[0], otherEnd ).index()
                    enOther = Edge3.edgeNumber[ver[0]][otherEnd]
                    verOther = tet.edgeMapping(enOther)
                    if verOther[0] == ver[0]:
                        adjacent = ( eiOther, seg )
                    else:
                        wtOther = surf.edgeWeight(eiOther).safeLongValue()
                        adjacent = ( eiOther, wtOther - seg )

                    # If the adjacent segment is one of the targets, then we
                    # are done; otherwise, we add it to the stack.
                    output = targets.get( adjacent, None )
                    if output is not None:
                        found.add(output)
                        #return output
                    else:
                        stack.append(adjacent)
            elif seg > f[0] + q:
                # The current segment belongs to a parallel triangular cell
                # at vertex ver[1].
                for otherEnd in range(4):
                    if otherEnd in { ver[0], ver[1] }:
                        continue

                    # The current segment is adjacent to a segment of the
                    # edge with endpoints ver[1] and otherEnd.
                    eiOther = tet.edge( ver[1], otherEnd ).index()
                    enOther = Edge3.edgeNumber[ver[1]][otherEnd]
                    verOther = tet.edgeMapping(enOther)
                    if verOther[0] == ver[1]:
                        adjacent = ( eiOther, wt - seg )
                    else:
                        wtOther = surf.edgeWeight(eiOther).safeLongValue()
                        adjacent = ( eiOther, wtOther - wt + seg )

                    # If the adjacent segment is one of the targets, then we
                    # are done; otherwise, we add it to the stack.
                    output = targets.get( adjacent, None )
                    if output is not None:
                        found.add(output)
                        #return output
                    else:
                        stack.append(adjacent)
            elif q > 0:
                # The quadrilaterals divide tet into two "sides". The edge
                # opposite this segment has endpoints lying on different
                # sides, so we label these endpoints accordingly.
                side = [ { 0, qType + 1 } ]
                side.append( {0,1,2,3} - side[0] )
                if ver[0] not in side[0]:
                    side[0], side[1] = side[1], side[0]
                side[0].remove( ver[0] )
                side[1].remove( ver[1] )
                opp = [ side[0].pop(), side[1].pop() ]

                # Find all edges containing segments that are adjacent to the
                # current segment.
                qDepth = seg - f[0]     # 0 <= qDepth <= q
                if qDepth == 0:
                    adjEndpoints = [ [ ver[0], opp[1] ] ]
                elif qDepth == q:
                    adjEndpoints = [ [ opp[0], ver[1] ] ]
                else:
                    adjEndpoints = [
                            [ ver[0], opp[1] ],
                            [ opp[0], ver[1] ],
                            [ opp[0], opp[1] ] ]
                for start, end in adjEndpoints:
                    # The current segment is adjacent to a segment of the
                    # edge going from start to end.
                    eiAdj = tet.edge( start, end ).index()
                    enAdj = Edge3.edgeNumber[start][end]
                    verAdj = tet.edgeMapping(enAdj)
                    triangles = surf.triangles(
                            teti, verAdj[0] ).safeLongValue()
                    if verAdj[0] == start:
                        adjacent = ( eiAdj, triangles + qDepth )
                    else:
                        adjacent = ( eiAdj, triangles + q - qDepth )

                    # If the adjacent segment is one of the targets, then we
                    # are done; otherwise, we add it to the stack.
                    output = targets.get( adjacent, None )
                    if output is not None:
                        found.add(output)
                        #return output
                    else:
                        stack.append(adjacent)

        # End of loop. Move on to the next embedding of edge e.

    # If the search terminates without finding the ideal edge, then the ideal
    # edge must belong to a component that gets destroyed.
    return found
    #return None


def _survivingSegments(surf):
    """
    Uses the given normal surface to divide the edges of the ambient
    triangulation into segments, and returns a dictionary describing the
    segments that would survive after crushing surf.

    In detail, the keys of the returned dictionary will be surviving segments
    encoded as pairs of the form (ei, s), where:
    --> ei is an edge index; and
    --> s is a segment number from 0 to w, inclusive, where w is the weight
        of surf on edge ei.
    The segments for each edge e are numbered in ascending order from the one
    incident to e.vertex(0) to the one incident to e.vertex(1). The returned
    dictionary will map each such segment to a pair (t, en), where:
    --> t is the index after crushing of a tetrahedron that will be incident
        to the segment in question; and
    --> en is an edge number (from 0 to 5, inclusive) of this tetrahedron
        that corresponds to the segment in question.
    """
    tri = surf.triangulation()
    survivors = dict()
    shift = 0
    for tet in tri.tetrahedra():
        teti = tet.index()
        #TODO Refactor using helper function.
        hasQuads = False
        for q in range(3):
            if surf.quads( teti, q ).safeLongValue() > 0:
                hasQuads = True
                break
        if hasQuads:
            # Presence of quads means tet is destroyed by crushing, which
            # will shift all larger tetrahedron indices down by one.
            shift += 1
            continue

        # No quads in tet, so there is a cell in the centre that survives
        # crushing. Find the edges of this cell that survive.
        for en in range(6):
            v = tet.edgeMapping(en)[0]
            ei = tet.edge(en).index()
            s = surf.triangles( teti, v ).safeLongValue()
            survivors[ (ei,s) ] = ( teti - shift, en )

    # Done!
    return survivors


if __name__ == "__main__":
    ## SFS [D: (2,1) (2,1)] U/m SFS [D: (2,1) (2,1)], m = [ 0,-1 | 1,0 ] : #1
    #sig = "hLALMkbcceffggemkbtibj"
    sig = argv[1]
    num = int( argv[2] )
    tri = Triangulation3.fromIsoSig(sig)
    qvsurfs = NormalSurfaces( tri, NormalCoords.NS_QUAD )
    surf = qvsurfs[num]
    #print(surf)

    #TODO
    survivors = _survivingSegments(surf)
    found = []
    for e in tri.edges():
        ei = e.index()
        wt = surf.edgeWeight(ei).safeLongValue()

        # Find ideal edges.
        for s in range( 1, wt ):
            result = findIdealEdges_old( surf, ( ei, s ), survivors )
            if result not in found:
                found.append(result)
    print()
    for f in found:
        print(f)

    #TODO Test _eulerWeights() routine.
    print()
    print( _eulerWeights(surf) )
