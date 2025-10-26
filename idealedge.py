"""
Find the ideal edges after crushing a normal surface.
"""
from regina import *
from loop import NotLoop, IdealLoop


def decomposeAlong( surf, oldLoops ):
    """
    Decomposes along surf, and returns a list of the resulting components.

    In detail, each item in the returned list is a list I of ideal loops,
    encoded as instances of IdealLoop, such that:
    --> each loop in I lies inside the same triangulation T; and
    --> the corresponding component is obtained by drilling out the loops in
        I from T.
    Thus, a side-effect of this routine is that it effectively deletes any
    components that contain no ideal loops at all (since there is no way to
    recover a triangulation from an empty list of ideal loops).

    Another side-effect is that this routine might detect and delete some
    ideal loops that are "trivial" in the sense that they bound embedded
    discs. However, note that this routine does not systematically test
    whether every loop is trivial or nontrivial, so it is still possible for
    the output to include some trivial loops.

    The given oldLoops list should be a list of pre-existing ideal loops,
    encoded as instances of IdealLoop.

    The given normal surface surf should be either:
    --> an annulus or 2-sphere that is disjoint from all of the pre-existing
        ideal loops; or
    --> a separating 2-sphere that intersects one of the pre-existing ideal
        loops in exactly two points, and is disjoint from all of the other
        pre-existing ideal loops.
    This routine raises ValueError if surf is not of one of these allowed
    types.

    We also require surf to be a quadrilateral vertex normal surface, but
    this routine does not check this condition.

    Pre-condition
    --> The given surf should be a quadrilateral vertex normal surface.
    --> If surf is an annulus, then each boundary component that it meets
        must be a two-triangle torus.
    --> The given surf and each ideal loop in oldLoops must all lie in the
        same triangulation.
    """
    #TODO WORKING HERE.
    # Find where the new ideal loops will be after crushing.
    #TODO Might need to update this once idealLoops() has been reimplemented.
    loopInfo = idealLoops( surf, oldLoops )
    crushed = surf.crush()

    # Split crushed into its components.
    if crushed.isConnected():
        components = [crushed]
        compLoopInfo = [loopInfo]
    else:
        components = list( crushed.triangulateComponents() )

        # Work out how tetrahedra get renumbered after splitting crushed into
        # its components.
        shiftedIndex = []
        compSize = [0] * crushed.countComponents()
        for i in range( crushed.size() ):
            compi = crushed.tetrahedron(i).component().index()
            shiftedIndex.append( compSize[compi] )
            compSize[compi] += 1

        # Shift tetrahedron indices for the ideal loops to account for the
        # renumbering computed above.
        compLoopInfo = [ [] for _ in range( crushed.countComponents() ) ]
        for seq in loopInfo:
            shiftedSeq = []
            #TODO This unpacking won't work anymore now that we're tracking
            #   orientations.
            for t, e in seq:
                compi = crushed.tetrahedron(t).component().index()
                shiftedSeq.append( ( shiftedIndex[t], e ) )
            compLoopInfo[compi].append(shiftedSeq)

    # Use compLoopInfo to find the ideal loops in each component.
    output = []
    for compi in range( crushed.countComponents() ):
        tri = components[compi]
        loopInfo = compLoopInfo[compi]
        loops = []
        for seq in loopInfo:
            edgeList = []
            #TODO This unpacking won't work anymore now that we're tracking
            #   orientations.
            for t, e in seq:
                edgeList.append( tri.tetrahedron(t).edge(e) )

            # Note that we could have a degenerate loop.
            try:
                loop = IdealLoop(edgeList)
            except NotLoop:
                # Ignore degenerate loop.
                continue
            else:
                loops.append(loop)
        if len(loops) == 1:
            try:
                loops[0].simplify()
                loops[0].simplify()
            except BoundsDisc:
                # Ignore trivial loop.
                continue
        output.append(loops)
    return output


def idealLoops( surf, oldLoops=[] ):
    """
    Returns information about the ideal loops after crushing the given normal
    surface surf.

    The given oldLoops list (which may be empty, and is empty by default)
    should be a list of pre-existing ideal loops, encoded as instances of
    IdealLoop. Each of these ideal loops must lie in the same triangulation
    as surf, and these ideal loops must all be mutually disjoint.

    The given normal surface surf should be either:
    --> an annulus or 2-sphere that is disjoint from all of the pre-existing
        ideal loops; or
    --> a separating 2-sphere that intersects one of the pre-existing ideal
        loops in exactly two points, and is disjoint from all of the other
        pre-existing ideal loops.
    This routine raises ValueError if surf is not of one of these allowed
    types.

    We also require surf to be a quadrilateral vertex normal surface, but
    this routine does not check this condition.

    This routine returns a list describing the ideal loops that would arise
    after crushing the given surface (see below for a more detailed
    description of how the ideal loops before crushing are related to the
    ideal loops after crushing). Each such ideal loop is encoded as a list of
    pairs of the form (ia, t, h), where:
    --> ia is the index after crushing of a tetrahedron that will be incident
        to one of the ideal edges;
    --> t is the vertex number (from 0 to 3, inclusive) of tetrahedron ia at
        the tail of the ideal edge in question; and
    --> h is the vertex number of tetrahedron ia at the head of the ideal
        edge.
    A caveat to this is that when the given surf is a 2-sphere, there is one
    possible degenerate ideal loop: a pair of edges giving an unknotted loop,
    such that the two edges get merged to become a single non-loop edge after
    crushing. This routine does not check for such degenerate loops, so they
    might appear in the returned list.

    Crushing the given surface has the following effects:
    --> Pre-existing ideal loops that are disjoint from the surface will be
        left topologically untouched. In particular, their orientations will
        be preserved.
    --> Ideal loops that intersect the surface will be split into multiple
        arcs, and each such arc may or may not survive to become a new ideal
        loop after crushing. The orientation will be preserved for the arcs
        that do survive.
    --> If the surface is an annulus (which, as specified above, must be
        disjoint from all pre-existing ideal loops), then crushing might
        create an entirely new ideal loop. This new loop will be assigned an
        arbitrary orientation.

    Pre-condition:
    --> The given surf should be a quadrilateral vertex normal surface.
    --> If surf is an annulus, then each boundary component that it meets
        must be a two-triangle torus.
    """
    # The given surf must be either a 2-sphere or an annulus. Moreover:
    # - In the 2-sphere case, we allow one of the ideal loops to have
    #   nonempty intersection with the surface.
    # - In the annulus case, we might create a new ideal loop by flattening
    #   a chain of boundary bigon faces.
    if isSphere(surf):
        loopMustBeDisjoint = False
        possibleLoopFromBoundary = False
    elif isAnnulus(surf):
        loopMustBeDisjoint = True
        possibleLoopFromBoundary = True
    else:
        allowed = "annuli and 2-spheres"
        msg = ( "This routine currently only accepts {} ".format(allowed) +
                "for the input surface." )
        raise ValueError(msg)

    # Find the ideal loops that arise from the pre-existing ideal loops.
    tri = surf.triangulation()
    newLoops = []
    targets = _survivingSegments(surf)
    for oldLoop in oldLoops:
        wt = oldLoop.weight(surf)
        if wt == 2:
            if loopMustBeDisjoint:
                msg = ( "Too many intersections between the surface and " +
                        "the pre-existing ideal loops." )
                raise ValueError(msg)

            # For a 2-sphere, we currently only allow at most one ideal loop
            # to intersect the surface.
            loopMustBeDisjoint = True
        elif wt != 0:
            msg = ( "Each ideal loop must intersect the surface in " +
                    "either exactly 0 points or exactly 2 points." )
            raise ValueError(msg)

        # The given surface splits the current oldLoop into some number of
        # components. Which of these components survive to become new ideal
        # loops after crushing?
        for arc in oldLoop.splitArcs(surf):
            seg = arc[0]
            idEdge = _findIdealEdge( surf, seg, targets )
            if idEdge is None:
                # This component does not survive after crushing.
                continue

            # This component survives after crushing.
            newLoop = [idEdge]
            for seg in arc[1:]:
                newLoop.append( _findIdealEdge( surf, seg, targets ) )
            newLoops.append(newLoop)

    # Will there also be an entirely new ideal loop created by flattening a
    # chain of boundary bigons?
    if possibleLoopFromBoundary and countIncidentBoundaries(surf) == 1:
        # Find a segment incident to the chain of boundary bigons.
        for e in tri.edges():
            ei = e.index()
            if ( e.isBoundary() and
                    surf.edgeWeight(ei).safeLongValue() >= 2 ):
                # Arbitrarily assign orientation +1.
                seg = ( ei, 1, 1 )
                break

        # If this segment survives after crushing, then it forms a new ideal
        # loop of length one.
        idEdge = _findIdealEdge( surf, seg, targets )
        if idEdge is not None:
            newLoops.append( [idEdge] )

    #TODO If we crushed an annulus, it would probably be useful to use
    #   fillIdealEdge() to include additional ideal loops obtained by filling
    #   in pinched 2-sphere boundary components.
    #
    #   If/when we implement this functionality, we will need to document the
    #   possibility that we could create an addiitonal new ideal loop. We
    #   should probably also note that this would come at the cost of
    #   introducing a new tetrahedron.

    # Done!
    return newLoops


def fillIdealEdge(tri):
    """
    If tri is an invalid component resulting from crushing an annulus, then
    closes this triangulation and returns the resulting ideal edge.

    This routine assumes that the boundary triangles of tri form a
    two-triangle 2-sphere with two of its three vertices pinched to form a
    single invalid vertex.
    """
    # Check some basic pre-conditions.
    if tri.isValid() or tri.countBoundaryTriangles() != 2:
        return

    # Might need to layer one tetrahedron to achieve the generic case.
    layerEdge = None
    for edge in tri.edges():
        if not edge.isBoundary():
            continue

        # Does this edge meet two distinct boundary triangles? If so, and if
        # it is also the *only* boundary edge with this property, then we
        # need to layer on this edge.
        emb = [ edge.embedding(0),
                edge.embedding( edge.degree() - 1 ) ]
        face = []
        for i in range(2):
            tet = emb[i].tetrahedron()
            faceNum = emb[i].vertices()[3-i]
            face.append( tet.triangle(faceNum) )
        if face[0] == face[1]:
            continue
        if layerEdge is None:
            layerEdge = edge
        else:
            layerEdge = None
            break
    if layerEdge is not None:
        tri.layerOn(layerEdge)

    # Find the ideal edge. It suffices to look for a boundary edge whose
    # endpoints are identified.
    for idEdge in tri.edges():
        if not idEdge.isBoundary():
            continue
        if idEdge.vertex(0) == idEdge.vertex(1):
            break

    # Finish by closing up the boundary.
    emb = [ idEdge.embedding(0),
            idEdge.embedding( idEdge.degree() - 1 ) ]
    ver = [ emb[i].vertices() for i in range(2) ]
    me = emb[0].tetrahedron()
    you = emb[1].tetrahedron()
    myFace = ver[0][3]
    gluing = Perm4(
            ver[0][0], ver[1][0],
            ver[0][1], ver[1][1],
            ver[0][2], ver[1][3],
            ver[0][3], ver[1][2] )
    me.join( myFace, you, gluing )

    # Now that the triangulation has changed, we can only access the ideal
    # edge through one of the tetrahedra in which it is embedded.
    return me.edge( ver[0][0], ver[0][1] )


def isAnnulus(s):
    """
    Is the given normal surface s an annulus?

    Pre-condition:
    --> It is known in advance that s is connected.
    """
    return ( s.isCompact() and s.isOrientable() and
            s.hasRealBoundary() and s.eulerChar() == 0 )


def isSphere(s):
    """
    Is the given normal surface s a 2-sphere?

    Pre-condition:
    --> It is known in advance that s is connected.
    """
    return ( s.isCompact() and s.isOrientable() and
            not s.hasRealBoundary() and s.eulerChar() == 2 )


def countIncidentBoundaries(s):
    """
    In the triangulation containing the given normal surface s, counts the
    number of boundary components that are incident to s.

    Pre-condition:
    --> The surface s lies inside a triangulation with only real boundary
        components.
    """
    tri = s.triangulation()
    incident = set()
    for e in tri.edges():
        bdy = e.boundaryComponent()
        if ( bdy is None ) or ( bdy.index() in incident ):
            continue
        if s.edgeWeight( e.index() ).safeLongValue() > 0:
            incident.add( bdy.index() )
    return len(incident)


def _findIdealEdge( surf, start, targets=None ):
    """
    Returns details of the ideal edge that corresponds to the given start
    segment after crushing surf.

    Specifically, if the ideal edge belongs to an ideal arc that gets
    destroyed after crushing, then this routine returns None. Otherwise, this
    routine returns a triple (i, t, h), where:
    --> i is the index after crushing of a tetrahedron that will be incident
        to the ideal edge;
    --> t is the vertex number of this tetrahedron that is at the tail of the
        ideal edge; and
    --> h is the vertex number that is at the head of the ideal edge.
    Here, "tail" and "head" are with respect to the orientation of the ideal
    edge, which will be consistent with the given start segment.

    If the dictionary of surviving segments has been precomputed using the
    _survivingSegments() routine, then this can be supplied using the
    optional targets argument. Otherwise, this routine will compute the
    surviving segments for itself.
    """
    tri = surf.triangulation()
    if targets is None:
        targets = _survivingSegments(surf)

    # If the start segment is one of the targets, then we are already done.
    output = targets.get( start, None )
    if output is not None:
        return output

    # Otherwise, we find the ideal edge using depth-first search.
    #
    # In theory, this could be done in polynomial time using the
    # Agol-Hass-Thurston weighted orbit-counting algorithm. However,
    # depth-first search is much easier to implement, and works very well in
    # practice.
    stack = [start]
    visited = set()
    while stack:
        current = stack.pop()
        if current in visited:
            continue

        # We haven't visited the current segment yet, so we need to find all
        # segments that are adjacent to it along parallel cells or faces.
        ei, seg, orient = current
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
                    # This is the quad type that is disjoint from e.
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
                    #
                    #           ver[0]
                    #              •
                    #             / \
                    #      edge e/   \
                    #           /     \
                    #    ver[1]•       •otherEnd
                    #
                    eiOther = tet.edge( ver[0], otherEnd ).index()
                    enOther = Edge3.edgeNumber[ver[0]][otherEnd]
                    verOther = tet.edgeMapping(enOther)
                    if verOther[0] == ver[0]:
                        # Same tails, hence same orientation.
                        adjacent = ( eiOther, seg, orient )
                    else:
                        # Opposite orientation.
                        wtOther = surf.edgeWeight(eiOther).safeLongValue()
                        adjacent = ( eiOther, wtOther - seg, -orient )

                    # If the adjacent segment is one of the targets, then we
                    # are done; otherwise, we add it to the stack.
                    output = targets.get( adjacent, None )
                    if output is not None:
                        return output
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
                    #
                    #           ver[0]
                    #              •
                    #             /
                    #      edge e/
                    #           /
                    #    ver[1]•-------•otherEnd
                    #
                    eiOther = tet.edge( ver[1], otherEnd ).index()
                    enOther = Edge3.edgeNumber[ver[1]][otherEnd]
                    verOther = tet.edgeMapping(enOther)
                    if verOther[0] == ver[1]:
                        # Opposite orientation.
                        adjacent = ( eiOther, wt - seg, -orient )
                    else:
                        # Same orientation.
                        wtOther = surf.edgeWeight(eiOther).safeLongValue()
                        adjacent = ( eiOther, wtOther - wt + seg, orient )

                    # If the adjacent segment is one of the targets, then we
                    # are done; otherwise, we add it to the stack.
                    output = targets.get( adjacent, None )
                    if output is not None:
                        return output
                    else:
                        stack.append(adjacent)
            elif q > 0:
                # At this point, we have f[0] <= seg <= f[0] + q.
                #
                # The quadrilaterals divide tet into two "sides". The edge
                # opposite this segment has endpoints lying on different
                # sides, and we can label these opposite endpoints opp[i],
                # i in {0,1}, so that ver[i] and opp[i] lie on the same side
                # of the quadrilaterals, as shown in the diagram below.
                #
                #               ver[0]
                #                  •
                #                 /|\
                #          edge e/ | \
                #               /__|__\
                #              /|  |  |\
                #       ver[1]•-|--|--|-•opp[1]
                #              \|__|__|/
                #               \  |  /
                #                \ | /opposite edge
                #                 \|/
                #                  •
                #               opp[0]
                #
                side = [ { 0, qType + 1 } ]
                side.append( {0,1,2,3} - side[0] )
                if ver[0] not in side[0]:
                    side[0], side[1] = side[1], side[0]
                side[0].remove( ver[0] )
                side[1].remove( ver[1] )
                opp = [ side[0].pop(), side[1].pop() ]

                # Find all edges containing segments that are adjacent to the
                # current segment.
                #
                # It is crucial that for each pair in adjEndpoints, the first
                # vertex is on the "0" side of the quadrilateral, and the
                # second vertex is on the "1" side of the quadrilateral.
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
                        # Same orientation.
                        adjacent = ( eiAdj, triangles + qDepth, orient )
                    else:
                        # Opposite orientation.
                        adjacent = ( eiAdj, triangles + q - qDepth, -orient )

                    # If the adjacent segment is one of the targets, then we
                    # are done; otherwise, we add it to the stack.
                    output = targets.get( adjacent, None )
                    if output is not None:
                        return output
                    else:
                        stack.append(adjacent)

        # End of loop. Move on to the next embedding of edge e.

    # If the search terminates without finding the ideal edge, then the ideal
    # edge must belong to a component that gets destroyed.
    return None


def _survivingSegments(surf):
    """
    Uses the given normal surface to divide the edges of the ambient
    triangulation into segments, and returns a dictionary describing the
    segments that would survive after crushing surf.

    In detail, the keys of the returned dictionary will be oriented surviving
    segments, encoded as triples of the form (ei, s, o), where:
    --> ei is an edge index;
    --> s is a segment number from 0 to w, inclusive, where w is the weight
        of surf on edge ei; and
    --> o is +1 if edge ei is oriented from vertex 0 to vertex 1, and -1 if
        edge ei is oriented from vertex 1 to vertex 0.
    The segments for each edge e are numbered in ascending order from the one
    incident to e.vertex(0) to the one incident to e.vertex(1). The returned
    dictionary will map each such segment to a triple (ia, t, h), where:
    --> ia is the index after crushing of a tetrahedron that will be incident
        to the segment in question;
    --> t is the vertex number (from 0 to 3, inclusive) of tetrahedron ia
        that is at the tail of the edge that corresponds to the segment in
        question; and
    --> h is the vertex number of tetrahedron ia that is at the head of the
        edge.
    """
    tri = surf.triangulation()
    survivors = dict()
    shift = 0
    for tet in tri.tetrahedra():
        teti = tet.index()
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
            tail = tet.edgeMapping(en)[0]
            head = tet.edgeMapping(en)[1]
            ei = tet.edge(en).index()
            s = surf.triangles( teti, tail ).safeLongValue()

            # Include both possible orientations.
            survivors[ (ei,s,1) ] = ( teti - shift, tail, head )
            survivors[ (ei,s,-1) ] = ( teti - shift, head, tail )

    # Done!
    return survivors
