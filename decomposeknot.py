"""
Decompose knots into prime knots.
"""
from regina import *
from idealedge import decomposeAlong


def embedInTriangulation(knot):
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
    if isinstance( knot, PacketOfLink ):
        packet = PacketOfTriangulation3(tri)
        packet.setLabel( knot.adornedLabel(
            "Embedded as edge #{}".format( tet.edge(edgeNum).index() ) ) )
        knot.insertChildLast(packet)
    return ( tri, tet.index(), edgeNum )
