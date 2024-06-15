"""
Uses an ideal loop given by 0/1 Dehn surgery to test whether a knot is
nontrivially knotted.
"""
from regina import *
from loop import IdealLoop


def surgery0(oldLoop):
    """
    Constructs a new ideal loop given by 0/1 Dehn surgery on the given ideal
    loop.

    Warning:
    --> This routine currently uses fast heuristics to attempt to construct
        the desired triangulation, and is not guaranteed to terminate.

    Returns:
        The newly constructed ideal loop.
    """
    # Triangulate the exterior with boundary edges appearing as the meridian
    # and longitude. The last step is not guaranteed to terminate in theory,
    # but it should be fine in practice.
    tri = oldLoop.drill()
    tri.idealToFinite()
    tri.intelligentSimplify()
    tri.intelligentSimplify()
    mer, lon = tri.meridianLongitude()

    # Get a tetrahedron index and edge number for the meridian, so that we
    # can remember its location after closing up the boundary.
    emb = mer.embedding(0)
    tet = emb.tetrahedron()
    edgeNum = emb.face()

    # Close up the boundary and build the new IdealLoop.
    layer = tri.layerOn(lon)
    layer.join( 0, layer, Perm4(0,1) )
    idealEdge = tet.edge(edgeNum)
    newLoop = IdealLoop( [idealEdge] )
    newLoop.simplify()
    return newLoop


def isKnotted(loop):
    """
    Is the given ideal loop nontrivially knotted?
    """
    core = surgery0(loop)
    while True:
        # INVARIANT:
        #   It is guaranteed that the original loop is unknotted if and only
        #   if there is a normal 2-sphere intersecting the current core loop
        #   in *exactly* one point.
        #
        # Our current task is to search for a quadrilateral vertex normal
        # 2-sphere that intersects the current core loop in *at most* one
        # point.
        # - If no such 2-sphere exists, then the original loop is
        #   nontrivially knotted.
        # - If we find a 2-sphere intersecting the current core loop in
        #   *exactly* one point, then the original loop is unknotted.
        # - If we find a 2-sphere disjoint from the current core loop, then
        #   we can crush to obtain a new core loop that still satisfies the
        #   invariant.
        tri = core.triangulation()
        enumeration = TreeEnumeration( tri, NS_QUAD )
        while True:
            # Get the next 2-sphere.
            #TODO
            pass
        #TODO
        pass
    #TODO
    return
