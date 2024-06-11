"""
Decompose knots into prime knots.
"""
from timeit import default_timer
from regina import *
from idealedge import decomposeAlong, isSphere
from loop import IdealLoop


def drill(loop):
    """
    Returns an ideal triangulation of the 3-manifold given by drilling out
    the given loop.
    """
    #TODO
    return


def embeddedLoopPacket(loop):
    """
    Returns a packet of the triangulation containing the given loop, with an
    ideal triangulation of the drilled 3-manifold as a child.
    """
    #TODO Update to use the drill() routine.
    tri = loop.triangulation()

    # Ideal triangulation given by drilling.
    drilled = PacketOfTriangulation3(tri)
    drillLocations = []
    for ei in loop:
        emb = tri.edge(ei).embedding(0)
        drillLocations.append( ( emb.tetrahedron(), emb.edge() ) )
    for tet, edgeNum in drillLocations:
        drilled.pinchEdge( tet.edge(edgeNum) )
    drilled.intelligentSimplify()
    drilled.setLabel( "Drilled: {}".format( drilled.isoSig() ) )

    # Edge-ideal triangulation, with drilled triangulation as child.
    packet = PacketOfTriangulation3(tri)
    packet.insertChildLast(drilled)
    return packet


def embedInTriangulation( knot, insertAsChild=False ):
    """
    Embeds the given knot as an ideal loop in a triangulation of the
    3-sphere.

    Warning:
    --> This routine currently uses fast heuristics to attempt to construct
        the desired triangulation, and is not guaranteed to terminate.
    """
    if knot.countComponents() > 1:
        raise ValueError( "Can only embed knots in a triangulation." )

    # Triangulate the exterior with boundary edges appearing as the meridian
    # and longitude. The last step is not guaranteed to terminate in theory,
    # but it should be fine in practice.
    tri = knot.complement()
    tri.intelligentSimplify()
    tri.idealToFinite()
    tri.intelligentSimplify()
    mer, lon = tri.meridianLongitude()

    # Get a tetrahedron index and edge number for the longitude, so that we
    # can remember its location after closing up the boundary.
    emb = lon.embedding(0)
    tet = emb.tetrahedron()
    edgeNum = emb.face()

    # Close up the boundary and build the IdealLoop.
    layer = tri.layerOn(mer)
    layer.join( 0, layer, Perm4(0,1) )
    idealEdge = tet.edge(edgeNum)
    loop = IdealLoop( [idealEdge] )
    if insertAsChild and isinstance( knot, PacketOfLink ):
        packet = embeddedLoopPacket(loop)
        packet.setLabel( knot.adornedLabel(
            "Embedded as edge {}".format( idealEdge.index() ) ) )
        knot.insertChildLast(packet)
    return loop


def decompose( knot, insertAsChild=False, timeout=10, verbose=False ):
    """
    Decomposes the given knot into prime pieces, represented as 3-spheres
    in which the prime knots are embedded as ideal loops.
    """
    start = default_timer()
    if verbose:
        prev = start
    primes = []
    if knot.complement().isSolidTorus():
        toProcess = []
    else:
        toProcess = [ embedInTriangulation(knot) ]
    steps = 0
    while toProcess:
        # INVARIANT:
        #   At this point, the following are guaranteed to hold:
        #   --> Each element of toProcess is an ideal loop forming a
        #       nontrivial knot.
        #   --> Each element of primes is an ideal loop forming a nontrivial
        #       *prime* knot.
        #   --> The input knot is given by composing all of the knots
        #       represented in toProcess and primes.
        oldLoop = toProcess.pop()
        tri = oldLoop.triangulation()

        # Search for a suitable quadrilateral vertex normal 2-sphere to
        # crush. If no such 2-sphere exists, then the oldLoop forms a
        # nontrivial prime knot.
        enumeration = TreeEnumeration( tri, NS_QUAD )
        while True:
            steps += 1
            time = default_timer()
            elapsed = time - start
            if elapsed > timeout:
                msg = "Timed out after {} step(s).".format(steps)
                raise RuntimeError(msg)
            elif verbose and time - prev > 5:
                prev = time
                msg = "Time: {:.6f}. Steps: {}.".format(
                        elapsed, steps )
                print(msg)

            # Get the next 2-sphere. If there are no more 2-spheres, then the
            # oldLoop must form a nontrivial prime knot.
            if enumeration.next():
                sphere = enumeration.buildSurface()
                if not isSphere(sphere):
                    continue
            else:
                primes.append(oldLoop)
                break

            #TODO Can minimise repeated code, and also tidy up the logic.
            # We only want 2-spheres that intersect the oldLoop in either
            # exactly 0 points or exactly 2 points.
            wt = oldLoop.weight(sphere)
            if wt == 0:
                # When the weight is 0, crushing does nothing topologically,
                # except possibly splitting off a single 3-sphere component
                # that doesn't even contain an ideal loop.
                decomposed = decomposeAlong( sphere, [oldLoop] )
                for newLoops in decomposed:
                    if newLoops:
                        toProcess.append( newLoops[0] )
                        break
                break
            elif wt != 2:
                continue

            # When the weight is 2, crushing does one of the following:
            # - The oldLoop gets decomposed into two new ideal loops, each of
            #   which forms a knot. Moreover, at least one of these knots
            #   must be nontrivial.
            # - We get a single new ideal loop that is topologically
            #   equivalent to the oldLoop, plus possibly a single extra
            #   3-sphere component that doesn't even contain an ideal loop.
            decomposed = decomposeAlong( sphere, [oldLoop] )
            knots = []
            for newLoops in decomposed:
                if newLoops:
                    knots.append( newLoop[0] )
            #TODO
            pass
        #TODO
        pass
    #TODO Account for multi-edge ideal loops.
            # At this point, we have a 2-sphere that intersects the ideal
            # edge in 2 points. Do we make progress after decomposing?
            decomposed = decomposeAlong( sphere, {edgeIndex} )
            nontrivialPieces = []
            nontrivialSize = 0
            for pieceTri, pieceLoops in decomposed:
                # Is the edge loop a nontrivial knot?
                pieceTetIndex, pieceEdgeNum = pieceLoops[0]
                drilled = Triangulation3(pieceTri)
                drilled.pinchEdge( drilled.tetrahedron(
                    pieceTetIndex ).edge(pieceEdgeNum) )
                drilled.intelligentSimplify()
                if drilled.isSolidTorus():
                    continue
                nontrivialPieces.append( (
                    pieceTri, pieceTetIndex, pieceEdgeNum ) )
                nontrivialSize += pieceTri.size()
            if ( nontrivialSize < tri.size()
                    or len(nontrivialPieces) > 1 ):
                # We made progress!
                toProcess.extend(nontrivialPieces)
                break
            else:
                #TODO This must be an error, otherwise algorithm doesn't
                #   actually work.
                raise RuntimeError( "Unexpected failure to make progress." )

    # Output some auxiliary information before returning the list of primes.
    elapsed = default_timer() - start
    msg = "Time: {:.6f}. Steps: {}.".format( elapsed, steps )
    print(msg)
    if insertAsChild and isinstance( knot, PacketOfLink ):
        container = Container( "Primes ({})".format(msg) )
        knot.insertChildLast(container)
        for i, p in enumerate(primes):
            tri, loop = p
            packet = embeddedLoopPacket( tri, loop )
            loopEdgeIndices = list(loop)
            if len(loop) == 1:
                adorn = "Embedded as edge {}".format( loopEdgeIndices[0] )
            else:
                indices = ""
                for ei in loopEdgeIndices[:-1]:
                    indices += ", {}".format(ei)
                adorn = "Embedded as edges {} and {}".format(
                        indices[2:], loopEdgeIndices[-1] )
            packet.setLabel( "Prime knot #{} ({})".format( i, adorn ) )
            container.insertChildLast(packet)
    return primes
