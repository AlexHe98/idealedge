"""
Decompose knots into prime knots.
"""
from regina import *
from idealedge import decomposeAlong


def embedInTriangulation(knot):
    """
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
