"""
Find the ideal edges after crushing a normal surface.
"""
from regina import *
from loop import IdealLoop  #TODO Probably need to keep this, but double-check.
from triloops import EdgeIdealTriangulation
from segment import OrientedSegment
from aux.bdry import hasOnlyMinimalRealTorusBoundaryComponents
from aux.tetrenum import tetRenumbering
from aux.quad import tetHasQuads
from aux.surface import SurfaceType, hasOnlyNonTrivialBoundaryCurves
from aux.surface import isSphere, isAnnulus, countIncidentBoundaries
from aux.looperror import NotLoop
from retriangulate.insert import layerOn


#TODO Overhaul decomposeAlong() to not only return the components that
#       survive crushing, but also to do the book-keeping of tracking deleted
#       components and counting the number of orbital compressions.
#TODO Update usage for new output format.
def decomposeAlong( surf, edgeIdealTri=None ):
    """
    Decomposes along surf, and returns a list of the resulting components.

    In detail, each item in the returned list will be instance of either
    EdgeIdealTriangulation or Triangulation3. The underlying triangulations
    will just be the components that result from crushing surf.triangulation()
    along surf. If there are any EdgeIdealTriangulation objects in the
    returned list, then the IdealLoop objects that are embedded in these
    triangulations indicate loops that should be drilled out to obtain
    topologically useful 3-manifolds obtained by decomposing along surf.

    If there are no pre-existing IdealLoop objects to track, then only surf
    needs to be supplied as input to this routine. Otherwise, both surf and
    edgeIdealTri should be supplied, in which case edgeIdealTri should be an
    instance of EdgeIdealTriangulation that tracks all of the pre-existing
    IdealLoop objects.

    This routine might (but is not guaranteed to) detect that some ideal loops
    are "trivial" in the sense that they bound embedded discs. This routine
    handles such a trivial ideal loop L in one of two ways:
    --> If L is the only ideal loop in its ambient triangulation, then this
        routine will silently delete L.
    --> Otherwise, if there are other ideal loops in the ambient triangulation
        of L, then this routine will raise BoundsDisc.

    A side-effect of this routine is that it might (but is not guaranteed to)
    detect and silently delete some (or all, or none) of the ideal loops that
    are "trivial" in the sense that they bound embedded discs.

    If there are no pre-existing ideal loops, then surf should be of one of
    the following types:
    --> A 2-sphere.
    --> A disc.
    --> An annulus with nontrivial boundary curves.
    --> A Mobius band with nontrivial boundary curve.
    Otherwise, letting W denote the weight of surf on the pre-existing ideal
    loops, surf should be of one of the following types:
    --> A 2-sphere with either W == 2 or W == 0.
    --> A disc with W == 1 and nontrivial boundary curve.
    --> A disc with W == 0.
    --> A projective plane with W == 1.
    This routine raises ValueError if surf is not of one of these allowed
    types.

    We also require surf to be a quadrilateral vertex normal surface, but
    this routine does not check this condition.

    Precondition
    --> The given surf should be a quadrilateral vertex normal surface.
    --> If surf has real boundary, then each boundary component that it meets
        must be a two-triangle torus.
    --> If edgeIdealTri is supplied, then edgeIdealTri.triangulation() should
        be the same as surf.triangulation(). In other words, edgeIdealTri and
        surf should both reference the same triangulation object in memory.
    """
    #TODO Decide whether to bother with allowing projective planes as input.

    #TODO Make a final decision on how to deal with trivial loops, and then
    #       check that trivial loops are actually dealt with as intended.

    # Find where the new ideal loops will be after crushing.
    #
    #NOTE The newIdealLoopEmbs() requires, and will check, that the given
    #       surface is of one of the allowed types.
    loopEmbSeqsInOldTri = newIdealLoopEmbs( surf, oldLoops )
    doomed = [ tet for tet in surf.triangulation().tetrahedra()
              if tetHasQuads( surf, tet.index() ) ]
    tetIndicesAfterCrush = tetRenumbering(doomed)
    crushed = surf.crush()
    loopEmbSeqs = []
    for oldEmbSequence in loopEmbSeqsInOldTri:
        embSequence = []
        for oldEmb in oldEmbSequence:
            crushedTet = crushed.tetrahedron(
                    tetIndicesAfterCrush[ oldEmb.tetrahedron().index() ] )
            embSequence.append( EdgeEmbedding3(
                crushedTet, oldEmb.vertices() ) )
        loopEmbSeqs.append(embSequence)

    # Split crushed into its components.
    if crushed.isConnected():
        components = [crushed]
        compLoopInfo = [[]]
        for embSequence in loopEmbSeqs:
            compLoopInfo[0].append(
                    [ ( emb.tetrahedron().index(), emb.vertices() )
                     for emb in embSequence ] )
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

        # Using the renumbering that we just computed, record shifted
        # tetrahedron indices for the ideal loops.
        compLoopInfo = [ [] for _ in range( crushed.countComponents() ) ]
        for embSequence in loopEmbSeqs:
            singleLoopInfo = []
            for survivingEmb in embSequence:
                teti = survivingEmb.tetrahedron().index()
                singleLoopInfo.append( ( shiftedIndex[teti],
                                        survivingEmb.vertices() ) )

            # Abuse the fact that teti persists beyond the scope of the above
            # for loop.
            compi = crushed.tetrahedron(teti).component().index()
            compLoopInfo[compi].append(singleLoopInfo)

    # Use compLoopInfo to find the ideal loops in each component.
    output = []
    for compi in range( crushed.countComponents() ):
        tri = components[compi]
        loopInfo = compLoopInfo[compi]
        loops = []
        for singleLoopInfo in loopInfo:
            # To construct an IdealLoop, we need:
            #   --> a list of edges, in order as we traverse the loop; and
            #   --> an orientation, which is either +1 if the first edge of
            #       the loop is oriented from vertex 0 to vertex 1, and -1 if
            #       the first edge is oriented from vertex 1 to vertex 0.
            edgeList = []
            for teti, ver in singleLoopInfo:
                edgeList.append(
                        tri.tetrahedron(teti).edge( ver[0], ver[1] ) )
            firstTet = tri.tetrahedron( singleLoopInfo[0][0] )
            firstVer = singleLoopInfo[0][1]
            firstEdgeNum = Edge3.faceNumber(firstVer)
            if firstVer[0] == firstTet.edgeMapping(firstEdgeNum)[0]:
                orientation = 1
            else:
                orientation = -1

            # Note that we could have a degenerate loop.
            try:
                loop = IdealLoop( edgeList, orientation )
            except NotLoop:
                #TODO Is ignoring these the right thing to do?

                # Ignore degenerate loop.
                continue
            else:
                loops.append(loop)

        # If we have any loops at all, then package them all together as a
        # single EdgeIdealTriangulation. Otherwise, just add an ordinary
        # Triangulation3 to the output list.
        #
        #TODO Consider simplifying before adding to the output list. But this
        #       requires deciding how this routine should behave if
        #       simplification raises BoundsDisc.
        if loops:
            output.append( EdgeIdealTriangulation(loops) )
        else:
            output.append(tri)

    # All done!
    return output


#TODO First draft of the documentation is based on the old newIdealLoopEmbs()
#       routine, and doesn't quite describe what this new routine will
#       eventually do, so we will need to update documentation accordingly.
def trackIdealSegments( surf, edgeIdealTri=None ):
    """
    Tracks ideal segments through the operation of crushing the given normal
    surface surf.

    This routine returns a pair of objects:
    (0) A list describing the ideal loops that would arise after crushing the
        given surface (see below for a more detailed description of how the
        ideal loops before crushing are related to the ideal loops after
        crushing). Each such ideal loop is encoded as a list of surviving
        edge embeddings.
    (1) An integer counting the number of orbital compression discs that
        would be cut along as a consequence of crushing the given surface.

    If there are no pre-existing IdealLoop objects to track, then only surf
    needs to be supplied as input to this routine. Otherwise, both surf and
    edgeIdealTri should be supplied, in which case edgeIdealTri should be an
    instance of EdgeIdealTriangulation that tracks all of the pre-existing
    IdealLoop objects.

    If there are no pre-existing ideal loops, then surf should be of one of
    the following types:
    --> A 2-sphere.
    --> A disc.
    --> An annulus with nontrivial boundary curves.
    --> A Mobius band with nontrivial boundary curve.
    Otherwise, letting W denote the weight of surf on the pre-existing ideal
    loops, surf should be of one of the following types:
    --> A 2-sphere with either W == 2 or W == 0.
    --> A disc with W == 1 and nontrivial boundary curve.
    --> A disc with W == 0.
    --> A projective plane with W == 1.
    Regardless of whether there are any pre-existing ideal loops, the ambient
    triangulation surf.triangulation() must either be closed or have minimal
    toroidal boundary. This routine raises ValueError if any of these
    conditions are not satisfied.

    We also require surf to be a quadrilateral vertex normal surface, but
    this routine does not check this condition.

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
    --> If edgeIdealTri is supplied, then edgeIdealTri.triangulation() should
        be the same as surf.triangulation(). In other words, edgeIdealTri and
        surf should both reference the same triangulation object in memory.
    """
    tri = surf.triangulation()
    if not hasOnlyMinimalRealTorusBoundaryComponents(tri):
        raise ValueError( "Triangulation is only allowed to have real " +
                         "boundary components that are two-triangle tori" )

    # Check that surf is of one of the required types.
    surfType = SurfaceType.recognise(surf)
    if edgeIdealTri is None:
        # No pre-existing ideal loops.
        if surfType == SurfaceType.SPHERE:
            #TODO
            raise NotImplementedError()
        elif surfType == SurfaceType.DISC:
            #TODO
            raise NotImplementedError()
        elif surfType == SurfaceType.ANNULUS:
            if not hasOnlyNonTrivialBoundaryCurves(surf):
                raise ValueError( "An annulus must have nontrivial " +
                                 "boundary curves" )
            #TODO
            raise NotImplementedError()
        elif surfType == SurfaceType.MOBIUS:
            if not hasOnlyNonTrivialBoundaryCurves(surf):
                raise ValueError( "A Mobius band must have nontrivial " +
                                 "boundary curve" )
            #TODO
            raise NotImplementedError()
        else:
            raise ValueError( "With no pre-existing ideal edges, we do " +
                             "not support {}".format(surfType) )
    else:
        # EdgeIdealTriangulation always holds a nonempty collection of ideal
        # loops.
        weight = edgeIdealTri.weight(surf)
        if surfType == SurfaceType.SPHERE:
            if weight not in {0, 2}:
                raise ValueError( "A 2-sphere must have ideal weight " +
                                 "either 0 or 2" )
            #TODO
            raise NotImplementedError()
        elif surfType == SurfaceType.DISC:
            if weight == 0:
                #TODO
                raise NotImplementedError()
            elif weight == 1:
                if not hasOnlyNonTrivialBoundaryCurves(surf):
                    raise ValueError( "A disc with ideal weight 1 must " +
                                     "have nontrivial boundary curve" )
                #TODO
                raise NotImplementedError()
            else:
                raise ValueError( "A disc must have ideal weight " +
                                 "either 0 or 1" )
        elif surfType == SurfaceType.RP3:
            if weight != 1:
                raise ValueError( "A projective plane must have ideal " +
                                 "weight 1" )
            #TODO
            raise NotImplementedError()
        else:
            raise ValueError( "With an edge-ideal triangulation, we do " +
                             "not support {}".format(surfType) )
    #TODO
    raise NotImplementedError()


#TODO Update documentation and implementation to:
#       --> use the new TriangulationWithEmbeddedLoops class, and
#       --> account for the extra cases that arise from SFS recognition.
def newIdealLoopEmbs( surf, oldLoops=[] ):
    """
    Returns surviving edge embeddings which describe the ideal loops after
    crushing the given normal surface surf.

    The given oldLoops list (which may be empty, and is empty by default)
    should be a list of pre-existing ideal loops, encoded as instances of
    IdealLoop. Each of these ideal loops must lie in the same triangulation
    as surf, and these ideal loops must all be mutually disjoint.

    If there are no pre-existing ideal loops, then surf should be of one of
    the following types:
    --> A 2-sphere.
    --> A disc.
    --> An annulus with nontrivial boundary curves.
    --> A Mobius band with nontrivial boundary curve.
    Otherwise, letting W denote the weight of surf on the pre-existing ideal
    loops, surf should be of one of the following types:
    --> A 2-sphere with either W == 2 or W == 0.
    --> A disc with W == 1 and nontrivial boundary curve.
    --> A disc with W == 0.
    --> A projective plane with W == 1.
    This routine raises ValueError if surf is not of one of these allowed
    types.

    We also require surf to be a quadrilateral vertex normal surface, but
    this routine does not check this condition.

    This routine returns a list describing the ideal loops that would arise
    after crushing the given surface (see below for a more detailed
    description of how the ideal loops before crushing are related to the
    ideal loops after crushing). Each such ideal loop is encoded as a list of
    surviving edge embeddings.

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
    #TODO Update the checks below to allow for the extra cases that arise in
    #       the context of bounded SFS recognition.

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
    newLoopEmbs = []
    survivors = OrientedSegment.survivors(surf)
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
        # arcs. Which of these arcs survive to become new ideal loops after
        # crushing?
        for arc in oldLoop.splitArcs(surf):
            seg = arc[0]
            survivingSeg = seg.translateAlongParallelCells(survivors)
            if survivingSeg is None:
                # This arc does not survive after crushing.
                continue

            # This arc survives after crushing.
            newLoop = [ survivingSeg.survivingEmbedding() ]
            for seg in arc[1:]:
                survivingSeg = seg.translateAlongParallelCells(survivors)
                newLoop.append( survivingSeg.survivingEmbedding() )
            newLoopEmbs.append(newLoop)

    # Will there also be an entirely new ideal loop created by flattening a
    # chain of boundary bigons?
    if possibleLoopFromBoundary and countIncidentBoundaries(surf) == 1:
        # Find a segment incident to the chain of boundary bigons.
        for e in tri.edges():
            ei = e.index()
            if ( e.isBoundary() and
                    surf.edgeWeight(ei).pythonValue() >= 2 ):
                # Arbitrarily assign orientation +1.
                seg = OrientedSegment( surf, ei, 1, 1 )
                break

        # If this segment survives after crushing, then it will form a new
        # ideal loop of length one.
        survivingSeg = seg.translateAlongParallelCells(survivors)
        if survivingSeg is not None:
            newLoopEmbs.append( [ survivingSeg.survivingEmbedding() ] )

    #TODO If we crushed an annulus, it would probably be useful to use
    #   fillIdealEdges() to include additional ideal loops obtained by filling
    #   in pinched 2-sphere boundary components.
    #
    #   If/when we implement this functionality, we will need to document the
    #   possibility that we could create an additional new ideal loop. We
    #   should probably also note that this would come at the cost of
    #   introducing a new tetrahedron.

    # Done!
    return newLoopEmbs


def fillIdealEdges( tri, endpoints ):
    """
    Fills in two-triangle boundary 2-spheres of tri that have a symmetric pair
    of vertices in the given set of endpoints, and returns a list of the
    resulting ideal edges.

    The endpoints set should be specified as a set of integers corresponding
    to vertex indices in the given triangulation tri.

    In detail, each two-triangle boundary 2-sphere B that we fill in falls
    under one of the following two cases:
    --> If B is isomorphic to the boundary of a triangular pillow, then any
        two out of the three vertices of B form a symmetric pair. Let e denote
        the edge of B that has both vertices in the given endpoints set. We
        fill in B by gluing the two faces of B together in the obvious way.
        The edge e becomes one of the ideal edges in the returned list.
    --> If B is isomorphic to the boundary of a snapped 3-ball, then the only
        symmetric pair of vertices is given by the two "poles" (that is, the
        two vertices not incident to the equator edge). We fill in B by
        attaching a snapped 3-ball. This introduces a new edge connecting the
        two poles, which becomes one of the ideal edges in the returned list.

    This routine modifies the given triangulation directly. If the
    triangulation is currently oriented, then the filling operation will
    preserve this orientation.
    """
    fillEdges = []  # Items will be triples ( tet, edgeNum, doLayer ).
    for bc in tri.boundaryComponents():
        # Is bc a real two-triangle 2-sphere boundary component?
        #
        # Note that bc might have vertices pinched together (tri is not
        # assumed to be valid), so we must build the Triangulation2 to be sure
        # of computing the Euler characteristic correctly.
        built = bc.build()
        if built.eulerChar() != 2:
            continue
        if not bc.isReal():
            continue
        if bc.size() != 2:
            continue

        # Do we have a pillow-boundary or a snapped-ball-boundary?
        equatorEdgeNum = None
        for e in range(3):
            if built.triangle(0).adjacentTriangle(e) == built.triangle(1):
                if equatorEdgeNum is None:
                    equatorEdgeNum = e
                else:
                    # For snapped-ball-boundary, the equator edge is the
                    # unique edge that the two boundary triangles have in
                    # common. Since we have now found two common edges, we
                    # cannot have snapped-ball-boundary.
                    equatorEdgeNum = None
                    break

        # Do we have a symmetric pair of vertices in the endpoints set?
        if equatorEdgeNum is None:
            # We have pillow-boundary, so the answer is "yes" if and only if
            # exactly two out of the three vertices are in the endpoints set.
            idealEdgeNum = None
            for e in range(3):
                # In bc.triangle(0), the ideal edge should be opposite the
                # unique vertex that is *not* in the endpoints set.
                if bc.triangle(0).vertex(e).index() not in endpoints:
                    if idealEdgeNum is None:
                        idealEdgeNum = e
                    else:
                        idealEdgeNum = None
                        break
            if idealEdgeNum is not None:
                # We already have the ideal edge, so doLayer should be False.
                emb = bc.triangle(0).edge(e).front()
                fillEdges.append(
                        ( emb.tetrahedron(), emb.face(), False ) )
        else:
            # We have snapped-ball-boundary, so the answer is "yes" if and
            # only if the "poles" (ie, the vertices not incident to the
            # equator edge) both lie in the endpoints set, and the vertex on
            # the equator edge does not lie in the endpoints set.
            if ( bc.triangle(0).edge(
                equatorEdgeNum ).vertex(0).index() in endpoints ):
                # Equator vertex is in the endpoints set.
                # Move on to the next boundary component.
                continue
            if bc.triangle(0).vertex(equatorEdgeNum).index() not in endpoints:
                # Pole vertex is not in the endpoints set.
                # Move on to the next boundary component.
                continue

            # We just need to check the pole incident to bc.triangle(1) now.
            otherPole = built.triangle(0).adjacentEdge(equatorEdgeNum)
            if bc.triangle(1).vertex(otherPole).index() in endpoints:
                # Later, we will need to layer across the equator edge to
                # obtain the ideal edge (ie, we need doLayer to be True).
                emb = bc.triangle(0).edge(equatorEdgeNum).front()
                fillEdges.append(
                        ( emb.tetrahedron(), emb.edge(), True ) )

    # Perform all the fillings.
    idealEdgeLocations = []
    for tet, edgeNum, doLayer in fillEdges:
        if doLayer:
            # Make sure to use oriented layering.
            tet = layerOn( tet.edge(edgeNum) )
            edgeNum = 5
        idealEdgeLocations.append( ( tet, edgeNum ) )
        idealEdge = tet.edge(edgeNum)

        # Close up the pillow boundary.
        #
        #       back    front
        #          0    0
        #          •    •
        #         /|    |\
        #        / |    | \
        #      3•  |    |  •2
        #        \ |    | /
        #         \|    |/
        #          •    •
        #          1    1
        #
        front = idealEdge.front()
        back = idealEdge.back()
        front.tetrahedron().join(
                front.vertices()[3],
                back.tetrahedron(),
                back.vertices() * Perm4(2,3) * front.vertices().inverse() )

    # All done!
    return [ tet.edge(edgeNum) for tet, edgeNum in idealEdgeLocations ]
