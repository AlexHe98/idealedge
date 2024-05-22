"""
Decompose knots into prime knots.
"""
from timeit import default_timer
from regina import *
from idealedge import decomposeAlong, isSphere


def embeddedLoopPacket( tri, tetIndex, edgeNum ):
    packet = PacketOfTriangulation3(tri)
    drilled = PacketOfTriangulation3(tri)
    packet.insertChildLast(drilled)
    drilled.pinchEdge( drilled.tetrahedron(tetIndex).edge(edgeNum) )
    drilled.intelligentSimplify()
    drilled.setLabel( "Drilled: {}".format( drilled.isoSig() ) )
    return packet


def embedInTriangulation( knot, insertAsChild=False ):
    """
    Constructs a triangulation of the 3-sphere in which the given knot is
    embedded as an edge loop.

    This routine returns a triple (t,i,e), where:
    --> t is a triangulation of the 3-sphere;
    --> i is the index of a tetrahedron incident to the knot edge; and
    --> e is an edge number of tetrahedron i that corresponds to the knot
        edge.

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

    # Close up the boundary.
    layer = tri.layerOn(mer)
    layer.join( 0, layer, Perm4(0,1) )
    if insertAsChild and isinstance( knot, PacketOfLink ):
        packet = embeddedLoopPacket( tri, tet.index(), edgeNum )
        packet.setLabel( knot.adornedLabel(
            "Embedded as edge #{}".format( tet.edge(edgeNum).index() ) ) )
        knot.insertChildLast(packet)
    return ( tri, tet.index(), edgeNum )


def decompose( knot, insertAsChild=False, timeout=10, verbose=False ):
    """
    Decomposes the given knot into prime pieces, represented as 3-spheres
    in which the prime knots are embedded as edge loops.
    """
    start = default_timer()
    if verbose:
        prev = start
    primes = []
    toProcess = [ embedInTriangulation(knot) ]
    steps = 0
    while toProcess:
        initial = toProcess.pop()
        tri, tetIndex, edgeNum = initial
        edgeIndex = tri.tetrahedron(tetIndex).edge(edgeNum).index()

        # Keep attempting to crush 2-spheres until we make progress, either
        # by reducing the total number of tetrahedra, or by decomposing into
        # two nontrivial knots. If this fails, then we have a prime piece.
        #NOTE In practice, enumerating in quad coordinates (without any extra
        #   constraints) seems to be much faster than enumerating in standard
        #   coordinates (both with and without a positive Euler constraint).
        enumeration = TreeEnumeration( tri, NS_QUAD )
        quadExhausted = False
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
            if not quadExhausted:
                #print("Quad search")
                if enumeration.next():
                    sphere = enumeration.buildSurface()
                    if not isSphere(sphere):
                        continue
                else:
                    quadExhausted = True
                    standard = TreeEnumeration_EulerPositive(
                            tri, NS_STANDARD )
            if quadExhausted:
                #print("Standard search")
                if standard.next():
                    sphere = standard.buildSurface()
                else:
                    primes.append(initial)
                    break

            # We only want 2-spheres that intersect the ideal edge in either
            # 0 points or 2 points.
            wt = sphere.edgeWeight(edgeIndex).safeLongValue()
            if wt == 0:
                # In this case, crushing does nothing topologically, except
                # possibly splitting off trivial 3-sphere components.
                decomposed = decomposeAlong( sphere, {edgeIndex} )
                for pieceTri, pieceLoops in decomposed:
                    if pieceLoops:
                        # Found the piece containing the knot.
                        break

                # Just to be sure, check that we made progress by reducing
                # the number of tetrahedra.
                if pieceTri.size() < tri.size():
                    pieceTetIndex, pieceEdgeNum = pieceLoops[0]
                    toProcess.append( (
                        pieceTri, pieceTetIndex, pieceEdgeNum ) )
                    break
                else:
                    continue
            if wt != 2:
                continue

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
    elapsed = default_timer() - start
    msg = "Time: {:.6f}. Steps: {}.".format( elapsed, steps )
    print(msg)
    if insertAsChild and isinstance( knot, PacketOfLink ):
        container = Container( "Primes ({})".format(msg) )
        knot.insertChildLast(container)
        for i, p in enumerate(primes):
            tri, tetIndex, edgeNum = p
            packet = embeddedLoopPacket( tri, tetIndex, edgeNum )
            packet.setLabel( "Prime knot #{} (Embedded as edge #{})".format(
                i, tri.tetrahedron(tetIndex).edge(edgeNum).index() ) )
            container.insertChildLast(packet)
    return primes
