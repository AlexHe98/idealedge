"""
Decompose knots into prime knots.
"""
from sys import stdout
from timeit import default_timer
from regina import *
from idealedge import decomposeAlong, isSphere
from loop import IdealLoop


def embeddedLoopPacket(loop):
    """
    Returns a packet of the triangulation containing the given loop, with an
    ideal triangulation of the drilled 3-manifold as a child.
    """
    drilled = PacketOfTriangulation3( loop.drill() )
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

    Returns:
        The constructed ideal loop.
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


def decompose( knot, verbose=False, insertAsChild=False ):
    """
    Decomposes the given knot into prime pieces, represented as 3-spheres
    in which the prime knots are embedded as ideal loops.
    """
    template = "        Time: {:.6f}. Steps: {}. Primes: {}. #Tri: {}."
    numTri = 0
    start = default_timer()
    if verbose:
        prev = start
        if insertAsChild:
            log = ""
    primes = []
    toProcess = [ embedInTriangulation(knot) ]
    steps = 0
    while toProcess:
        # INVARIANT:
        #   At this point, the following are guaranteed to hold:
        #   --> Each element of toProcess is an ideal loop forming a knot.
        #   --> Each element of primes is an ideal loop forming a nontrivial
        #       prime knot.
        #   --> The input knot is given by composing all of the knots
        #       represented in toProcess and primes.
        oldLoop = toProcess.pop()
        tri = oldLoop.triangulation()
        numTri += 1
        if verbose:
            prev = default_timer()
            msg = "Edge-ideal: "
            if tri.size() == 1:
                msg += "1 tetrahedron.\n"
            else:
                msg += "{} tetrahedra.\n".format( tri.size() )
            msg += template.format(
                    prev - start, steps, len(primes), numTri )
            print(msg)
            stdout.flush()
            if insertAsChild:
                log += msg + "\n"

        # Search for a suitable quadrilateral vertex normal 2-sphere to
        # crush. If no such 2-sphere exists, then the oldLoop is prime.
        enumeration = TreeEnumeration( tri, NS_QUAD )
        while True:
            steps += 1
            time = default_timer()
            if verbose and time - prev > 5:
                prev = time
                msg = template.format(
                        time - start, steps, len(primes), numTri )
                print(msg)
                stdout.flush()
                if insertAsChild:
                    log += msg + "\n"

            # Get the next 2-sphere.
            if enumeration.next():
                sphere = enumeration.buildSurface()
                if not isSphere(sphere):
                    continue
            else:
                # No suitable 2-sphere means oldLoop is prime. But we only
                # care about the case where this prime is nontrivial.
                if verbose:
                    # Let the user know that we are about to check
                    # nontriviality, which is typically a bottleneck.
                    time = default_timer()
                    msg = template.format(
                            time - start, steps, len(primes), numTri )
                    msg += "\nFound prime knot! Is it nontrivial?"
                    drilled = oldLoop.drill()
                    drilled.idealToFinite()
                    drilled.intelligentSimplify()
                    drilled.intelligentSimplify()
                    msg += "\nDrilled: "
                    if drilled.size() == 1:
                        msg += "1 tetrahedron."
                    else:
                        msg += "{} tetrahedra.".format( drilled.size() )
                    print(msg)
                    stdout.flush()
                    if insertAsChild:
                        log += msg + "\n"

                    # Perform the nontriviality check.
                    if drilled.isSolidTorus():
                        ans = "\nNo, it's the unknot."
                    else:
                        ans = "\nYes, found a nontrivial prime knot!"
                        primes.append(oldLoop)

                    # Let the user know that we made it out the other side.
                    prev = default_timer()
                    msg = template.format(
                            prev - start, steps, len(primes), numTri )
                    msg += ans
                    print(msg)
                    stdout.flush()
                    if insertAsChild:
                        log += msg + "\n"
                elif not oldLoop.drill().isSolidTorus():
                    primes.append(oldLoop)
                break

            # We only want 2-spheres that intersect the oldLoop in either
            # exactly 0 points or exactly 2 points, since crushing such a
            # 2-sphere either:
            # - simplifies the triangulation containing the ideal loop;
            # - decomposes the oldLoop into two simpler knots; or
            # - (if oldLoop is unknotted) destroys all traces of the loop.
            wt = oldLoop.weight(sphere)
            if wt != 0 and wt != 2:
                continue
            decomposed = decomposeAlong( sphere, [oldLoop] )
            knots = []
            for newLoops in decomposed:
                if newLoops:
                    # We are guaranteed to have len(newLoops) == 1.
                    knots.append( newLoops[0] )
            for newLoop in knots:
                toProcess.append(newLoop)
            if verbose:
                prev = default_timer()
                msg = template.format(
                        prev - start, steps, len(primes), numTri )
                print(msg)
                stdout.flush()
                if insertAsChild:
                    log += msg + "\n"
            break

    # Output some auxiliary information before returning the list of primes.
    msg = template.format(
            default_timer() - start, steps, len(primes), numTri )
    print(msg)
    stdout.flush()
    if insertAsChild and isinstance( knot, PacketOfLink ):
        if verbose:
            log += msg
            container = Text(log)
        else:
            container = Container()
        container.setLabel( "Primes ({})".format( msg.lstrip() ) )
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
