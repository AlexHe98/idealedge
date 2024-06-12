"""
Decompose knots into prime knots.
"""
from sys import stdout
from timeit import default_timer
from threading import Timer
from regina import *
from idealedge import decomposeAlong, isSphere
from loop import IdealLoop


def drill(loop):
    """
    Returns an ideal triangulation of the 3-manifold given by drilling out
    the given loop.
    """
    drilled = Triangulation3( loop.triangulation() )
    drillLocations = []
    for ei in loop:
        emb = drilled.edge(ei).embedding(0)
        drillLocations.append( ( emb.tetrahedron(), emb.edge() ) )
    for tet, edgeNum in drillLocations:
        drilled.pinchEdge( tet.edge(edgeNum) )
    drilled.intelligentSimplify()
    return drilled


def embeddedLoopPacket(loop):
    """
    Returns a packet of the triangulation containing the given loop, with an
    ideal triangulation of the drilled 3-manifold as a child.
    """
    drilled = PacketOfTriangulation3( drill(loop) )
    drilled.setLabel( "Drilled: {}".format( drilled.isoSig() ) )
    packet = PacketOfTriangulation3( loop.triangulation() )
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
    tri.intelligentSimplify()
    #TODO
    #mer, lon = tri.meridianLongitude()
    def terminator( signum, frame ):
        raise TimeoutError()
    timeout = 1
    temp = Triangulation3(tri)
    while True:
        print( "    Attempt {}.".format(timeout) )
        stdout.flush()
        alarm = Timer( timeout, terminator )
        try:
            temp.meridianLongitude()
        except TimeoutError:
            timeout += 1
            temp = Triangulation3(tri)
            temp.subdivide()
            temp.intelligentSimplify()
            temp.intelligentSimplify()
        else:
            alarm.cancel()
            break
    tri = temp
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
    loop.simplify()
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
        #TODO
        if verbose:
            print( "    {}-vertex, {}-tet.".format(
                tri.countVertices(), tri.size() ) )
            stdout.flush()

        # Search for a suitable quadrilateral vertex normal 2-sphere to
        # crush. If no such 2-sphere exists, then the oldLoop forms a
        # nontrivial prime knot.
        enumeration = TreeEnumeration( tri, NS_QUAD )
        while True:
            steps += 1
            time = default_timer()
            elapsed = time - start
            if elapsed > timeout:
                msg = "Timed out ({:.6f}) after {} step(s).".format(
                        elapsed, steps )
                raise RuntimeError(msg)
            elif verbose and time - prev > 5:
                prev = time
                msg = "Time: {:.6f}. Steps: {}. Primes: {}.".format(
                        elapsed, steps, len(primes) )
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

            # We only want 2-spheres that intersect the oldLoop in either
            # exactly 0 points or exactly 2 points, since crushing such a
            # 2-sphere either:
            # - simplifies the triangulation containing the ideal loop; or
            # - decomposes the oldLoop into two simpler knots.
            wt = oldLoop.weight(sphere)
            if wt != 0 and wt != 2:
                continue
            decomposed = decomposeAlong( sphere, [oldLoop] )
            knots = []
            for newLoops in decomposed:
                if newLoops:
                    knots.append( newLoops[0] )
            if len(knots) == 1:
                toProcess.append( knots[0] )
                break

            # In the case where we decomposed oldLoop into simpler knots, we
            # only want to keep the nontrivial pieces.
            for newLoop in knots:
                if not drill(newLoop).isSolidTorus():
                    toProcess.append(newLoop)
            break

    # Output some auxiliary information before returning the list of primes.
    elapsed = default_timer() - start
    msg = "Time: {:.6f}. Steps: {}. Primes: {}.".format(
            elapsed, steps, len(primes) )
    print(msg)
    if insertAsChild and isinstance( knot, PacketOfLink ):
        container = Container( "Primes ({})".format(msg) )
        knot.insertChildLast(container)
        for i, primeLoop in enumerate(primes):
            packet = embeddedLoopPacket(primeLoop)
            loopEdgeIndices = list(primeLoop)
            if len(primeLoop) == 1:
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
