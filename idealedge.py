"""
Find the ideal edge after crushing a separating annulus in a one-vertex
triangulation of a 3-manifold with torus boundary.
"""
from timeit import default_timer
from regina import *
from loop import NotLoop, IdealLoop


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


def survivingSegments(surf):
    """
    Uses the given normal surface to divide the edges of the ambient
    triangulation into segments, and returns a dictionary describing the
    segments that would survive after crushing surf.
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


def _findIdealEdge( surf, start, targets=None ):
    """
    Returns details of the ideal edge that corresponds to the given start
    segment after crushing surf.

    Specifically, if the ideal edge belongs to a component that gets
    destroyed after crushing, then this routine returns None. Otherwise, this
    routine returns a pair (t,e), where t is the index (after crushing) of a
    tetrahedron that will meet the ideal edge after crushing, and e is an
    edge number of this tetrahedron that corresponds to the ideal edge.
    """
    tri = surf.triangulation()
    if targets is None:
        targets = survivingSegments(surf)

    # If the start segment is one of the targets, then we are already done.
    output = targets.get( start, None )
    if output is not None:
        return output

    # Otherwise, we find ideal edge using depth-first search.
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
                        return output
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
                        return output
                    else:
                        stack.append(adjacent)

        # End of loop. Move on to the next embedding of edge e.

    # If the search terminates without finding the ideal edge, then the ideal
    # edge must belong to a component that gets destroyed.
    return None


def idealLoops( surf, oldLoops=[] ):
    """
    Returns information about the ideal loops after crushing surf.

    The given oldLoops list (which is empty by default) should be a list of
    pre-existing ideal loops, encoded as instances of IdealLoop.

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
    after crushing the given surface. Each such ideal loop is encoded as a
    list of pairs of the form (t,e), where:
    --> t is the index (after crushing) of a tetrahedron that will meet one
        of the ideal edges; and
    --> e is an edge number of this tetrahedron that corresponds to the ideal
        edge in question.
    A caveat to this is that when the given surf is a 2-sphere, there is one
    possible degenerate ideal loop: a pair of edges giving an unknotted loop,
    such that the two edges get merged to become a single non-loop edge after
    crushing. This routine does not check for such degenerate loops, so they
    might appear in the returned list.

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
    targets = survivingSegments(surf)
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

        # What happens to the current oldLoop after crushing?
        for comp in oldLoop.components(surf):
            seg = comp[0]
            idEdge = _findIdealEdge( surf, seg, targets )
            if idEdge is None:
                # This component does not survive after crushing.
                continue

            # This component survives after crushing.
            newLoop = [idEdge]
            for seg in comp[1:]:
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
                seg = ( ei, 1 )
                break
        idEdge = _findIdealEdge( surf, seg )
        if idEdge is not None:
            newLoops.append( [idEdge] )

    # Done!
    return newLoops


def printAnnulusIdealEdges(surfaces):
    """
    Prints details of all ideal edges obtained by crushing annuli in the
    given list of vertex normal surfaces.

    Pre-condition:
    --> For each annulus in the given list of surfaces, every boundary
        component incident to this annulus must be a two-triangle torus.
    """
    for i, surf in enumerate(surfaces):
        if isAnnulus(surf):
            print( i, idealLoops(surf) )


def printSphereIdealEdges( surfaces, idealEdgeIndex ):
    """
    Prints details of all ideal edges obtained by crushing 2-spheres in the
    given list of vertex normal surfaces.
    """
    for i, surf in enumerate(surfaces):
        if isSphere(surf):
            if surf.edgeWeight(idealEdgeIndex).safeLongValue() != 2:
                continue
            idealLoop = IdealLoop( [
                surf.triangulation().edge(idealEdgeIndex) ] )
            print( i, idealLoops( surf, [idealLoop] ) )


def fillIdealEdge(tri):
    """
    If tri is an invalid component resulting from crushing an annulus, then
    closes this triangulation and returns the ideal edge.

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


def recogniseSummands( tri, threshold=40 ):
    """
    Attempts to recognise the prime summands of the given triangulation.

    This routine only proceeds with performing the prime decomposition if the
    number of tetrahedra in tri is strictly less than the threshold (default
    40), and returns True if and only if this is the case.
    """
    if tri.size() >= threshold:
        return False
    summands = tri.summands()
    if len(summands) == 0:
        tri.setLabel( tri.label() + ": S3" )
    elif len(summands) == 1:
        # Try combinatorial recognition.
        std = StandardTriangulation.recognise( summands[0] )
        if std is None:
            name = "Prime, not recognised"
        else:
            name = std.manifold().name()
        tri.setLabel( tri.label() + ": {}".format(name) )
    else:
        tri.setLabel( tri.label() + ": Non-prime" )

        # Find *all* quad vertex normal 2-spheres.
        surfs = NormalSurfaces( tri, NS_QUAD, NS_VERTEX )
        sphereFilter = SurfaceFilterProperties()
        sphereFilter.setEulerChars( [2] )
        sphereFilter.setCompactness( BoolSet(True) )
        sphereFilter.setOrientability( BoolSet(True) )
        sphereFilter.setRealBoundary( BoolSet(False) )
        spheres = PacketOfNormalSurfaces( surfs, sphereFilter )
        spheres.setLabel( "Quad vertex 2-spheres (Total: {})".format(
            spheres.size() ) )
        tri.insertChildLast(spheres)

        # Classify the summands.
        sumContainer = Container( "Summands (Total: {})".format(
            len(summands) ) )
        tri.insertChildLast(sumContainer)
        for sumNum, s in enumerate(summands):
            summand = PacketOfTriangulation3(s)
            sumContainer.insertChildLast(summand)

            # Try to combinatorially recognise this summand.
            std = StandardTriangulation.recognise(summand)
            if std is None:
                name = "Prime, not recognised"
            else:
                name = std.manifold().name()
            summand.setLabel( "Summand #{}: {}".format(
                sumNum, name ) )
    return True


def snapEdge(edge):
    """
    If the endpoints of the given edge are distinct, then uses a snapped ball
    to pinch these two endpoints together.

    This operation is equivalent to performing the following two operations:
    (1) Pinching the edge, which introduces a two-tetrahedron with a single
        degree-one edge e at its heart.
    (2) Performing a 2-1 edge move on e.

    If the triangulation containing the given edge is currently oriented,
    then this operation will preserve the orientation.

    Pre-conditions:
    --> The given edge belongs to a triangulation with no boundary faces.

    TODO:
    --> The stated pre-conditions are stronger than necessary.

    Parameters:
    --> edge    The edge whose endpoints should be snapped together.

    Returns:
        True if and only if snapping the given edge is possible.
    """
    if edge.vertex(0) == edge.vertex(1):
        return False
    tri = edge.triangulation()
    tri.pinchEdge(edge)

    # To find the degree-one edge at the heart of the pinch edge gadget, look
    # at the last two tetrahedra in tri.
    found = False
    for tetIndex in [ tri.size() - 1, tri.size() - 2 ]:
        for edgeNum in range(6):
            e = tri.tetrahedron(tetIndex).edge(edgeNum)
            if e.degree() == 1:
                found = True
                break
        if found:
            break

    # Finish up by performing a 2-1 move on e.
    if not tri.twoOneMove( e, 0 ):
        if not tri.twoOneMove( e, 1 ):
            raise RuntimeError( "Snap edge failed unexpectedly." )
    return True


def decomposeAlong( surf, oldLoops ):
    """
    Decomposes along surf, and returns a list of the resulting components.

    In detail, each item in the returned list is a list I of ideal loops,
    encoded as instances of IdealLoop, such that:
    --> each loop in I lies inside the same triangulation T; and
    --> the corresponding component is obtained by drilling out the loops in
        I from T.

    The given oldLoops list should be a nonempty list of pre-existing ideal
    loops, encoded as instances of IdealLoop.

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
    loopInfo = idealLoops( surf, oldLoops )
    crushed = surf.crush()

    # Split crushed into its components.
    if crushed.isConnected():
        components = [crushed]
        compLoopInfo = [loopInfo]
    else:
        components = list( crushed.triangulateComponents() )

        # Work out what the new tetrahedron indices are in all the components
        # of crushed.
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
            for t, e in seq:
                edgeList.append( tri.tetrahedron(t).edge(e) )

            # Note that we could have a degenerate loop.
            try:
                loop = IdealLoop(edgeList)
            except NotLoop:
                continue
            else:
                loops.append(loop)
        #TODO Simplify if possible.
        output.append(loops)
    return output


def decomposeAlongSpheres( surfaces, idealEdgeIndex, threshold=30 ):
    """
    """
    results = Container( "Decompose along 2-spheres" )
    surfaces.insertChildLast(results)
    for surfNum, surf in enumerate(surfaces):
        if not isSphere(surf):
            continue
        try:
            #TODO This needs to be updated.
            pieces = decomposeAlong( surf, {idealEdgeIndex} )
        except ValueError:
            continue
        container = Container( "Decompose along #{}".format(surfNum) )
        results.insertChildLast(container)
        for i, piece in enumerate(pieces):
            tri = PacketOfTriangulation3( piece[0] )
            loops = piece[1]
            tri.setLabel( "Component #{}: {}".format(
                i, loops ) )
            container.insertChildLast(tri)

            # Is tri a 3-sphere?
            if ( tri.knowsSphere() or tri.size() < threshold ):
                if tri.isSphere():
                    name = "S3"
                else:
                    name = "Not S3"
            else:
                name = "Not recognised"
            tri.setLabel( tri.label() + ": {}".format(name) )

            # Build drilled 3-manifold.
            drilled = PacketOfTriangulation3(tri)
            drilled.setLabel( tri.adornedLabel(
                "Pinched ideal edges" ) )
            tri.insertChildLast(drilled)
            for t, e in loops:
                drilled.pinchEdge( drilled.tetrahedron(t).edge(e) )
                drilled.intelligentSimplify()

                # Is drilled a solid torus?
                if ( drilled.knowsSolidTorus() or
                        drilled.size() < threshold ):
                    if drilled.isSolidTorus():
                        name = "Ideal solid torus"
                    else:
                        name = "Ideal, not solid torus"
                else:
                    name = "Ideal, not recognised"
                drilled.setLabel(
                        drilled.label() + ": {}".format(name) )
            #TODO
            pass
        #TODO
        pass
    #TODO
    return


def crushAnnuli( surfaces, threshold=30 ):
    """
    Crushes all annuli in the given packet of normal surfaces, and adds a
    Container of the resulting triangulations as a child of the given packet.

    This routine attempts to identify the topology of the manifold that
    results from crushing. The main strategy is to simplify and attempt
    combinatorial recognition. Additionally, whenever this routine encounters
    a component whose number of tetrahedra is strictly less than the
    threshold (default 30), it will also use more computationally intensive
    recognition algorithms involving normal surfaces.

    Pre-condition:
    --> For each annulus in the given list of surfaces, every boundary
        component incident to this annulus must be a two-triangle torus.
    """
    start = default_timer()
    results = Container( "Crush annuli" )
    surfaces.insertChildLast(results)
    annulusCount = 0
    for surfNum, surf in enumerate(surfaces):
        if not isAnnulus(surf):
            continue
        annulusCount += 1
        print()
        print( "Time: {:.6f}. Crush #{}.".format(
            default_timer() - start, surfNum) )

        # Is the current annulus a thin edge link?
        thin = surf.isThinEdgeLink()
        if thin[0] is None:
            thinAdorn = ""

        # Crush, and find the ideal edge amongst the components of the
        # resulting triangulation.
        tri = PacketOfTriangulation3( surf.crush() )
        tri.setLabel( "Crushed #{}".format(surfNum) )
        thin = surf.isThinEdgeLink()
        if thin[0] is not None:
            # Adorn label with details of this thin edge link.
            adorn = "Thin edge {}".format( thin[0].index() )
            if thin [1] is not None:
                adorn += " and {}".format( thin[1].index() )
            tri.setLabel( tri.adornedLabel(adorn) )
        components = []
        results.insertChildLast(tri)
        idEdgeDetails = idealLoops(surf)
        if idEdgeDetails:
            # There is only one ideal loop, given by a length-1 sequence of
            # ideal edges.
            idEdge = idEdgeDetails[0][0]
        else:
            idEdge = None
        idComp = None
        if tri.isEmpty():
            tri.setLabel( tri.label() + ": Empty" )
        else:
            if tri.isConnected():
                components.append(tri)
                if idEdge is not None:
                    idComp = 0
            else:
                tri.setLabel( tri.label() + ": Disconnected" )
                for compNum, c in enumerate( tri.triangulateComponents() ):
                    comp = PacketOfTriangulation3(c)
                    comp.setLabel( "Component #{}".format(compNum) )
                    components.append(comp)
                    tri.insertChildLast(comp)

                # Find the component containing the ideal edge, and adjust
                # the ideal tetrahedron index.
                if idEdge is not None:
                    idComp = tri.tetrahedron( idEdge[0] ).component().index()
                    idTeti = 0
                    for tet in tri.tetrahedra():
                        if tet.component().index() == idComp:
                            if tet.index() == idEdge[0]:
                                idEdge = ( idTeti, idEdge[1] )
                                break
                            else:
                                idTeti += 1

        # Go through the components and try to identify their topology.
        for compNum, comp in enumerate(components):
            print( "    Time: {:.6f}. Component #{}.".format(
                default_timer() - start, compNum ) )
            if not comp.isValid():
                comp.setLabel( comp.label() + ": INVALID" )

                # Fill in invalid boundary.
                filled = PacketOfTriangulation3(comp)
                invIdEdge = fillIdealEdge(filled)
                filled.setLabel( comp.adornedLabel(
                    "Closed, ideal edge {}".format(
                        invIdEdge.index() ) ) )
                comp.insertChildLast(filled)

                # Just in case, let's see if we can simplify and identify the
                # manifold given by drilling out the ideal edge.
                drilled = PacketOfTriangulation3(filled)
                filled.insertChildLast(drilled)
                ide = drilled.edge( invIdEdge.index() )
                drilled.setLabel( comp.adornedLabel(
                    "Closed, pinched edge {}".format( ide.index() ) ) )
                drilled.pinchEdge(ide)
                drilled.intelligentSimplify()
                if ( ( drilled.knowsSolidTorus() or
                    drilled.size() < threshold ) and
                    drilled.isSolidTorus() ):
                    name = "Ideal solid torus"
                else:
                    # Try to combinatorially recognise after truncating the
                    # ideal vertex.
                    trunc = PacketOfTriangulation3(drilled)
                    drilled.insertChildLast(trunc)
                    trunc.idealToFinite()
                    trunc.intelligentSimplify()
                    std = StandardTriangulation.recognise(trunc)
                    if std is None:
                        name = "Not recognised"
                        if drilled.knowsSolidTorus():
                            name += ", not solid torus"
                    else:
                        name = std.manifold().name()
                    trunc.setLabel( drilled.adornedLabel(
                        "Truncated" ) + ": {}".format(name) )
                drilled.setLabel(
                        drilled.label() + ": {}".format(name) )

                # Decompose the filled manifold into prime pieces (unless it
                # has too many tetrahedra).
                print( "        Attempted prime decomposition: {}.".format(
                    recogniseSummands( filled, threshold ) ) )
            else:
                # If this component contains the ideal edge, then attempt to
                # simplify (and possibly identify) the drilled manifold.
                if compNum == idComp:
                    drilled = PacketOfTriangulation3(comp)
                    comp.insertChildLast(drilled)
                    ide = drilled.tetrahedron( idEdge[0] ).edge( idEdge[1] )

                    # Need to label *before* drilling.
                    drilled.setLabel( comp.adornedLabel(
                        "Pinched edge {}".format( ide.index() ) ) )
                    comp.setLabel( comp.adornedLabel(
                        "Ideal edge {}".format( ide.index() ) ) )
                    drilled.pinchEdge(ide)
                    drilled.intelligentSimplify()

                    # Try to recognise the drilled manifold.
                    if ( ( drilled.knowsSolidTorus() or
                        drilled.size() < threshold ) and
                        drilled.isSolidTorus() ):
                        name = "Ideal solid torus"
                    else:
                        # Try to combinatorially recognise after truncating
                        # the ideal vertex.
                        trunc = PacketOfTriangulation3(drilled)
                        drilled.insertChildLast(trunc)
                        trunc.idealToFinite()
                        trunc.intelligentSimplify()
                        std = StandardTriangulation.recognise(trunc)
                        if std is None:
                            name = "Not recognised"
                            if drilled.knowsSolidTorus():
                                name += ", not solid torus"
                        else:
                            name = std.manifold().name()
                        trunc.setLabel( drilled.adornedLabel(
                            "Truncated" ) + ": {}".format(name) )
                    drilled.setLabel(
                            drilled.label() + ": {}".format(name) )

                # Decompose this component into prime pieces (unless this
                # component has too many tetrahedra).
                print( "        Attempted prime decomposition: {}.".format(
                    recogniseSummands( comp, threshold ) ) )

    # All done!
    print()
    print( "Time: {:.6f}. All done!".format(
        default_timer() - start ) )
    results.setLabel( results.adornedLabel(
        "Total {}".format(annulusCount) ) )
