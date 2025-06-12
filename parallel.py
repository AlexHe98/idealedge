"""
Find the ideal edges after crushing a normal surface.
"""
from sys import argv
from regina import *
from loop import NotLoop, IdealLoop


_quadSameSide = [
        Perm4(1,0,3,2),
        Perm4(2,3,0,1),
        Perm4(3,2,1,0) ]
_quadOpposite = [
        Perm4(2,3,0,1),
        Perm4(1,0,3,2),
        Perm4(1,0,3,2)]


def parallelityBaseTopology(surf):
    """
    Determines the topology of the base B of each component of the
    parallelity bundle by computing the Euler characteristic and the number
    of boundary curves of B.
    """
    weights = _eulerWeights(surf)
    parOrbits = _computeParallelityOrbits(surf)
    survivingSegments = _survivingSegments(surf)
    parBdries = _parallelityBoundaries( surf, survivingSegments )

    # The _eulerWeights() routine assigns segment weights so that the total
    # weight is *twice* the Euler characteristic of the base.
    euler = [ wt//2 for wt in
             _parallelityBundleWeights( surf, weights, parOrbits ) ]

    # Work out which parallelity boundaries belong to each orbit.
    bdryCurves = []
    for orbit in parOrbits:
        bdryCurves.append( parBdries.intersection(orbit) )

    # Put information together.
    output = []
    for i, eulerChar in enumerate(euler):
        boundaries = bdryCurves[i]
        survivors = dict()
        for b in boundaries:
            survivors[b] = survivingSegments.get( b, None )
        output.append( ( eulerChar, survivors ) )
    return output


def _parallelityBoundaries( surf, survivingSegments=None ):
    """
    Finds all boundary components of the parallelity bundle.

    In detail, each such boundary component B is identified using a single
    segment S inside B; in the case where B contains at least one surviving
    segment, the chosen segment S is guaranteed to be one of the surviving
    segments. This routine returns a set consisting of these chosen segments,
    with exactly one such segment for each boundary component of the
    parallelity bundle.
    """
    tri = surf.triangulation()
    if survivingSegments is None:
        survivingSegments = _survivingSegments(surf)
    parallelBoundaries = set()

    # Find where type-1 and type-2 segments change between thick and thin
    # regions, since we will need this information to be able to traverse
    # boundary components of the parallelity bundle.
    regionChanges = []
    for edge in tri.edges():
        regionChanges.append( _segmentRegionChanges( edge, surf ) )

        # Special case
        # ------------
        # We will consider a component of the parallelity bundle consisting
        # entirely of a single segment S to have boundary given by S itself;
        # in other words, S is an "isolated" parallelity segment. Such
        # isolated parallelity segments can be characterised as central
        # segments that have no changes between thick and thin regions.
        isolatedPos = _centralSegment( edge.embedding(0), surf )
        if ( isolatedPos is None ) or ( isolatedPos in regionChanges[-1] ):
            continue
        parallelBoundaries.add( ( edge.index(), isolatedPos ) )

    # Traverse vertical boundary components of the parallelity boundary until
    # we have visited every boundary parallel face.
    parBdryFaceSegEmbs = _parallelityBoundaryFaceSegmentEmbeddings(surf)
    while parBdryFaceSegEmbs:
        _, startSegEmbs = parBdryFaceSegEmbs.popitem()
        currentSegEmb, stopSegEmb = startSegEmbs
        representativeSeg = ( currentSegEmb[0], currentSegEmb[2] )
        isRepSegSurviving = ( representativeSeg in survivingSegments )

        # Traverse the current boundary component.
        while True:
            currentEdgeInd, currentEmbInd, currentSegPos = currentSegEmb

            # Walk around the link of the current segment until we find the
            # next boundary parallel face.
            currentChanges = regionChanges[currentEdgeInd][currentSegPos]
            for changeInd, changeData in enumerate(currentChanges):
                change, embInd, _ = changeData
                if embInd == currentEmbInd:
                    break
            if change > 0:
                # We are currently at a change from thin to thick, so we need
                # to walk around in the forwards direction.
                changeInd += 1
                if changeInd == len(currentChanges):
                    changeInd = 0
            else:
                # We are currently at a change from thick to thin, so we need
                # to walk around in the backwards direction.
                changeInd -= 1
            nextRegionChange = currentChanges[changeInd]
            _, nextEmbInd, nextParFace = nextRegionChange
            currentSegEmb = ( currentEdgeInd, nextEmbInd, currentSegPos )

            # Have we come back to the first parallel face in this boundary
            # component of the parallelity bundle?
            if currentSegEmb == stopSegEmb:
                break

            # Continue traversing along this boundary component of the
            # parallelity bundle by walking along the next parallel face.
            nextParFaceSegEmbs = parBdryFaceSegEmbs.pop(nextParFace)
            i = nextParFaceSegEmbs.index(currentSegEmb)
            currentSegEmb = nextParFaceSegEmbs[1-i]

            # If necessary, update the representative segment.
            if isRepSegSurviving:
                continue
            currentSeg = ( currentSegEmb[0], currentSegEmb[2] )
            if currentSeg in survivingSegments:
                representativeSeg = currentSeg
                isRepSegSurviving = True

        # Record the representative segment for the current boundary
        # component of the parallelity bundle.
        parallelBoundaries.add(representativeSeg)

    # All done!
    return parallelBoundaries


def _centralSegment( edgeEmb, surf ):
    """
    Returns the position of the type-1 or type-2 central segment for the
    given edge embedding (with respect to the given normal surface), or None
    if there is no such segment.
    """
    # To have a central segment, the tetrahedron must not contain a quad
    # intersecting the edge.
    teti = edgeEmb.tetrahedron().index()
    eNum = edgeEmb.face()
    quadType = tetQuadType( surf, teti )
    if ( quadType is not None ) and ( quadType not in {eNum,5-eNum} ):
        return None

    # Central segment occurs immediately after all the triangles at vertex 0
    # of the edge.
    # We don't want to include type-0 segments, so check that there is at
    # least one triangle at either vertex 0 or vertex 1 of the edge.
    verts = edgeEmb.vertices()
    triCount = [ surf.triangles( teti, verts[i] ).safeLongValue()
                for i in range(2) ]
    if triCount[0] > 0 or triCount[1] > 0:
        return triCount[0]
    return None


def _segmentRegionChanges( edge, surf ):
    """
    Records the places where type-1 and type-2 segments of the given edge
    (with respect to the given normal surface) change regions from either
    thick to thin or thin to thick.

    In detail, this routine returns a dictionary that maps each type-1 or
    type-2 segment S to a list L that specifies the region changes for S. The
    entries of such a list L are of the form ( c, i, (t,f) ), where:
    --> c is +1 if the change is from thin to thick, and -1 if the change is
        from thick to thin;
    --> i is the index of the edge embedding at which the change occurs; and
    --> (t,f) specifies the parallel face that witnesses the region change,
        using a tetrahedron index t and a face number f.
    If the segment has no changes between thick and thin regions, then the
    segment will not appear as a key in the returned dictionary.
    """
    edgeEmbIndices = _edgeEmbeddingIndices( surf.triangulation() )
    edgeWt = surf.edgeWeight( edge.index() ).safeLongValue()
    regionChanges = dict()
    for embIndex, emb in enumerate( edge.embeddings() ):
        edgeNum = emb.face()
        ver = emb.vertices()
        tet = emb.tetrahedron()
        teti = tet.index()

        # Changes from thick to thin or thin to thick occur when tet contains
        # a quad that intersects the given edge.
        quadType = tetQuadType( surf, teti )
        if ( quadType is None ) or ( quadType in {edgeNum,5-edgeNum} ):
            continue

        # Changes occur at the segments identified below.
        changeSegs = []
        triCount = [ surf.triangles( teti, ver[i] ).safeLongValue()
                    for i in range(2) ]
        changeSegs.append( triCount[0] )
        changeSegs.append( edgeWt - triCount[1] )

        # The following permutations ensure that for any value of v in
        # {0,1,2,3}, the vertex labels satisfy the following:
        #
        #                  sameSide[v]
        #                       *
        #                      /|\
        #                     / | \
        #                    /  |  \
        #                   /   |   \
        #                  /____|____\
        #                 /|    |    |\
        #                / |    |    | \
        #   opposite[v] *--|----|----|--* opposite[sameSide[v]]
        #                \ |    |    | /
        #                 \|___/|\___|/
        #                  \  / | \  /
        #                   \/__|__\/
        #                    \  |  /
        #                     \ | /
        #                      \|/
        #                       *
        #                       v
        #
        sameSide = _quadSameSide[quadType]
        #opposite = _quadOpposite[quadType]

        # Record the changes.
        for j in range(2):
            seg = changeSegs[j]
            if seg is None:
                continue
            v = ver[j]

            # Change is from thick to thin if sameSide[v] == ver[2], and from
            # thin to thick if sameSide[v] == ver[3].
            if sameSide[v] == ver[2]:
                change = -1
            else:
                change = 1

            # Record the index of the edge embedding where the change occurs,
            # together with the location of the parallel face that witnesses
            # this change.
            parFace = ( teti, sameSide[v] )
            data = ( change, embIndex, parFace )
            if seg in regionChanges:
                regionChanges[seg].append(data)
            else:
                regionChanges[seg] = [data]

    # All done!
    return regionChanges


def _edgeEmbeddingIndices(tri):
    embIndices = { teti: { n: None for n in range(6) }
                  for teti in range( tri.size() ) }
    for edge in tri.edges():
        ei = edge.index()
        for embIndex, emb in enumerate( edge.embeddings() ):
            teti = emb.tetrahedron().index()
            edgeNum = emb.face()
            embIndices[teti][edgeNum] = ( ei, embIndex )
    return embIndices


def _parallelityBoundaryFaceSegmentEmbeddings(surf):
    """
    Returns information about how the two vertical edges of each parallelity
    boundary face are embedded as segments with respect to the given normal
    surface.

    In detail, each parallelity boundary face P is specified using a pair of
    the form (t,f), where:
    --> t is the index of the tetrahedron incident to P on the "outside"; and
    --> f is the face number of tetrahedron t specifying the triangular face
        that contains P.
    The output of this routine is a dictionary that maps each such pair to a
    2-element list L, where each element of L specifies an embedding of one
    of one of the vertical edges of P as a segment. Specifically, each such
    segment embedding is given by a tuple of the form (e,i,s), where:
    --> e is the index of the edge containing the segment;
    --> i is the index of the edge embedding; and
    --> s is the position of the segment along edge e.
    """
    tri = surf.triangulation()
    edgeEmbeddingIndices = _edgeEmbeddingIndices(tri)
    faceSegments = dict()
    for tet in tri.tetrahedra():
        teti = tet.index()
        quadType = tetQuadType( surf, teti )
        if quadType is None:
            continue

        # Find parallelity boundary faces incident to tet.
        sameSide = _quadSameSide[quadType]
        opposite = _quadOpposite[quadType]
        # The above permutations ensure that for any value of face in
        # {0,1,2,3}, the vertex labels satisfy the following:
        #
        #                        face
        #                          *
        #                         /|\
        #                        / | \
        #                       /  |  \
        #                      /   |   \
        #                     /____|____\
        #                    /|    |    |\
        #                   / |    |    | \
        #   opposite[face] *--|----|----|--* opposite[sameSide[face]]
        #                   \ |    |    | /
        #                    \|___/|\___|/
        #                     \  / | \  /
        #                      \/__|__\/
        #                       \  |  /
        #                        \ | /
        #                         \|/
        #                          *
        #                   sameSide[face]
        for face in range(4):
            # We have a corner or parallel face P that is boundary (or
            # isolated, in the case where it is "boundary on both sides").
            parFace = ( teti, face )
            faceSegments[parFace] = []
            endpts = [ face, sameSide[face] ]
            triCount = surf.triangles( teti, sameSide[face] ).safeLongValue()
            for i in range(2):
                # Look at one of the edges incident to P.
                edgeNum = Edge3.faceNumber( Perm4(
                        sameSide[face], opposite[endpts[i]],
                        face, opposite[endpts[1-i]] ) )
                edgeIndex = tet.edge(edgeNum).index()
                embIndex = edgeEmbeddingIndices[teti][edgeNum][1]

                # Find the segment incident to P.
                if tet.edgeMapping(edgeNum)[0] == sameSide[face]:
                    segPosition = triCount
                else:
                    edgeWt = surf.edgeWeight(edgeIndex).safeLongValue()
                    segPosition = edgeWt - triCount
                faceSegments[parFace].append(
                        ( edgeIndex, embIndex, segPosition ) )

    return faceSegments


def _parallelityBundleWeights( surf, weights, parOrbits ):
    """
    Calculates the total weight of each component of the parallelity bundle.

    In detail, this routine returns a list of orbit weights, in the same
    order as the orbits in the given parOrbits list.

    WARNING:
    --> This routine does not check whether the given segment weights make
        sense for the given normal surface. Therefore, this routine may
        produce undefined behaviour when given invalid input arguments.
    """
    orbitWeights = []
    for orbit in parOrbits:
        totalWeight = 0
        for seg in orbit:
            totalWeight += weights[seg]
        orbitWeights.append(totalWeight)
    return orbitWeights


def _computeParallelityOrbits(surf):
    # Naively compute orbits using repeated depth-first search.
    # At least in theory, it would be better to do this using something like
    # the Agol-Hass-Thurston weighted orbit-counting algorithm.
    segments = _segments(surf)
    parOrbits = []
    while segments:
        seg, segType = segments.popitem()
        if segType == 0:
            continue

        # Compute the orbit of this type-1 or type-2 segment.
        orbit = _computeOrbit( surf, seg )
        for otherSeg in orbit:
            if otherSeg != seg:
                del segments[otherSeg]
        parOrbits.append(orbit)
    return parOrbits


def _computeOrbit( surf, start ):
    tri = surf.triangulation()

    # Depth-first search.
    stack = [start]
    orbit = set()
    while stack:
        current = stack.pop()
        if current in orbit:
            continue
        orbit.add(current)

        # We haven't visited the current segment yet, so we need to find all
        # segments that are adjacent to it along corner/parallel cells or
        # corner/parallel faces.
        ei, seg = current
        e = tri.edge(ei)
        wt = surf.edgeWeight(ei).safeLongValue()
        for emb in e.embeddings():
            tet = emb.tetrahedron()
            teti = tet.index()
            en = emb.face()
            ver = emb.vertices()

            # To locate the relevant corner/parallel cells and
            # corner/parallel faces in tet, need to get the normal
            # coordinates incident to e.
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

            # Does the current segment belong to a corner/parallel cell or
            # corner/parallel face in this tet? If so, then we need to find
            # all adjacent segments.
            if seg < f[0]:
                # The current segment belongs to a corner/parallel triangular
                # cell at vertex ver[0].
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
                    stack.append(adjacent)
            elif seg > f[0] + q:
                # The current segment belongs to a corner/parallel triangular
                # cell at vertex ver[1].
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
                    stack.append(adjacent)

        # End of loop. Move on to the next embedding of edge e.

    # All done!
    return orbit


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

    # Vertices of the base correspond to type-1 or type-2 segments.
    weights = dict()
    segments = _segments(surf)
    for seg in segments:
        segType = segments[seg]
        if segType != 0:
            weights[seg] = 2
        else:
            # Could set weights for type-0 segments, but this is not
            # necessary.
            #weights[seg] = 0
            pass

    # Edges of the base correspond to corner/parallel faces of the induced
    # cell decomposition, and faces of the base correspond to corner/parallel
    # cells.
    #
    # Let:  E = number of edges of the base
    #       F = number of faces of the base
    #       q = number of parallel quad cells
    #       t = number of corner/parallel triangle cells
    #       b = number of corner/parallel faces in the vertical boundary
    #       i = number of isolated corner/parallel faces 
    # We adjust weights using the following observation:
    #   -2E + 2F = -(4q + 3t + b + 2i) + 2(q + t) = -2q - t - b - 2i.
    for tet in tri.tetrahedra():
        teti = tet.index()

        # For each corner/parallel triangle cell, adjust the weight on one
        # incident segment.
        for triType in range(4):
            triCount = surf.triangles( teti, triType ).safeLongValue()
            if triCount == 0:
                continue

            # We have at least one corner/parallel triangle cell.
            if triType in {0,3}:
                incidentEdgeNum = 2
            else:
                incidentEdgeNum = 3
            ei = tet.edge(incidentEdgeNum).index()
            if tet.edgeMapping(incidentEdgeNum)[0] == triType:
                for segment in range( 0, triCount ):
                    weights[ (ei,segment) ] -= 1
            else:
                edgeWt = surf.edgeWeight(ei).safeLongValue()
                for segment in range( edgeWt, edgeWt - triCount, -1 ):
                    weights[ (ei,segment) ] -= 1

        # We also need to adjust weights using parallel quad cells, and
        # parallel faces that are either boundary or isolated. Such cases
        # only arise in tetrahedra that contain at least one quad.
        quads = tetQuads( surf, teti )
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

        # For each corner/parallel face that is either isolated or on the
        # boundary of the parallelity bundle, adjust the weight on one
        # incident segment. In the isolated case, we only need to subtract 1
        # on this side, because we will eventually subtract another 1 on the
        # other side as well.
        sameSide = _quadSameSide[quadType]
        opposite = _quadOpposite[quadType]
        for face in range(4):
            # A corner/parallel face that is isolated/boundary occurs
            # immediately after the triangles (if any) at the relevant corner
            # of the tetrahedron.
            triCount = surf.triangles( teti, sameSide[face] ).safeLongValue()
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

    # All done!
    return weights


def tetQuadType( surf, tetIndex ):
    """
    Returns the quad type in which the given normal surface intersects the
    tetrahedron with the given index, or None if there is no such quad.
    """
    quads = tetQuads( surf, tetIndex )
    if quads is None:
        return None
    else:
        return quads[0]


def tetQuads( surf, tetIndex ):
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
    dictionary will map each surviving segment to a pair (t, en), where:
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
        if tetQuadType( surf, teti ) is not None:
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

    # Test parallelityBaseTopology() routine.
    print()
    print( "parallelityBaseTopology(surf)" )
    parBaseTop = parallelityBaseTopology(surf)
    for b in parBaseTop:
        print(b)

    # Print some extra information that is helpful for manual experimentation
    # in the GUI.
    print()
    print( "Ideal edges" )
    crushed = surf.crush()
    comp = crushed.triangulateComponents()
    compSize = [ 0 for _ in range( crushed.countComponents() ) ]
    compTeti = []
    for tet in crushed.tetrahedra():
        compi = tet.component().index()
        compTeti.append( compSize[compi] )
        compSize[compi] += 1
    for eulerChar, parBdryData in parBaseTop:
        if eulerChar == 1:
            # Don't need to worry about contractible parallelity components,
            # and we can assume that no parallelity components have
            # projective plane base.
            continue

        # Print details of ideal edges.
        for parBdrySeg, survivor in parBdryData.items():
            if survivor is None:
                continue

            # Surviving edge is an ideal edge.
            idTeti, idEdgeNum = survivor
            idTet = crushed.tetrahedron(idTeti)
            idCompi = idTet.component().index()
            idCompTeti = compTeti[idTeti]
            idEdge = comp[idCompi].tetrahedron(idCompTeti).edge(idEdgeNum)
            print( "{}: Component #{}, Edge #{}".format(
                parBdrySeg, 1 + idCompi, idEdge.index() ) )
