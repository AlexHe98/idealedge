"""
Uses an ideal loop given by 0/1 Dehn surgery to test whether a knot is
nontrivially knotted.
"""
from regina import *
from loop import IdealLoop


def surgery0(knot):
    """
    Constructs an ideal loop given by 0/1 Dehn surgery on the given knot.

    Warning:
    --> This routine currently uses fast heuristics to attempt to construct
        the desired triangulation, and is not guaranteed to terminate.

    Returns:
        The constructed ideal loop.
    """
    if knot.countComponents() > 1:
        raise ValueError( "Can only perform the surgery on a knot." )

    # Triangulate the exterior with boundary edges appearing as the meridian
    # and longitude. The last step is not guaranteed to terminate in theory,
    # but it should be fine in practice.
    tri = knot.complement()
    tri.intelligentSimplify()
    tri.idealToFinite()
    tri.intelligentSimplify()
    tri.intelligentSimplify()
    mer, lon = tri.meridianLongitude()

    # Get a tetrahedron index and edge number for the meridian, so that we
    # can remember its location after closing up the boundary.
    emb = mer.embedding(0)
    tet = emb.tetrahedron()
    edgeNum = emb.face()

    # Close up the boundary and build the IdealLoop.
    layer = tri.layerOn(lon)
    layer.join( 0, layer, Perm4(0,1) )
    idealEdge = tet.edge(edgeNum)
    loop = IdealLoop( [idealEdge] )
    loop.simplify()
    return loop


def isKnotted(knot):
    """
    Is the given knot nontrivially knotted?
    """
    #TODO
    return
