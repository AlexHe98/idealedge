"""
Uses an ideal loop given by 0/1 Dehn surgery to test whether a knot is
nontrivially knotted.
"""
from regina import *
from idealedge import decomposeAlong, isSphere
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


def isKnotted( loop, tracker=None ):
    """
    Is the given ideal loop nontrivially knotted?

    The tracker is an optional DecompositionTracker.
    """
    core = surgery0(loop)
    while True:
        # INVARIANT:
        #   It is guaranteed that the original loop is unknotted if and only
        #   if there is a normal 2-sphere intersecting the current core loop
        #   in *exactly* one point.
        tri = core.triangulation()
        if tracker is not None:
            tracker.newTri( tri.size() )

        # Let L denote the current core loop. If there is a normal 2-sphere
        # intersecting L in exactly one point, then we have the following
        # sequence of implications:
        #   ==> There is a standard vertex normal 2-sphere meeting L in
        #       either exactly zero points or exactly one point.
        #   ==> After adding positive Euler characteristic as an extra linear
        #       constraint, there is still a vertex normal 2-sphere meeting L
        #       in either exactly zero points or exactly one point.
        #
        # With this in mind, we now proceed by searching for vertex surfaces
        # with the positive Euler characteristic constraint included.
        # - If every such surface intersects L in two or more points, then
        #   the original loop is nontrivially knotted.
        # - If we find a projective plane that intersects L in exactly one
        #   point (assuming the original loop forms a knot in the 3-sphere,
        #   note that it is not possible to find a projective plane disjoint
        #   from L), then the original loop is nontrivially knotted.
        # - If we find a 2-sphere that intersects L in exactly one point,
        #   then the original loop is unknotted.
        # - If we find a 2-sphere that is disjoint from L, then we can try to
        #   crush to obtain a new core loop that still satisfies the
        #   invariant. This might fail because it destroys an edge in the
        #   core loop, but such a failure would certify that the original
        #   loop is unknotted.
        enumeration = TreeEnumeration_EulerPositive( tri, NS_STANDARD )
        while True:
            if tracker is not None:
                tracker.newSearch()

            # Get the next surface.
            if enumeration.next():
                surface = enumeration.buildSurface()
                wt = core.weight(surface)
                if wt > 1:
                    continue
                elif not isSphere(surface):
                    # Projective plane with weight one implies that the
                    # original loop is nontrivially knotted.
                    return True
            else:
                # As above, no suitable surface means the original loop is
                # nontrivially knotted.
                return True

            # At this point, we have a 2-sphere with weight at most one.
            if wt == 1:
                # The original loop is unknotted.
                return False
            else:
                # As above, when the weight is zero, we try crushing. If
                # anything goes wrong, then we certify that the original loop
                # is unknotted.
                try:
                    decomposed = decomposeAlong( surface, [core] )
                except Exception as e:
                    print( "Certified unknotted by exception: {}".format(e) )
                    return False
                else:
                    foundNewLoop = False
                    for newLoops in decomposed:
                        if not newLoops:
                            continue
                        elif len(newLoops) == 1:
                            if len( newLoops[0] ) != len(core):
                                print( "Certified unknotted by " +
                                        "shortened loop." )
                                return False
                            core = newLoops[0]
                            foundNewLoop = True
                            break
                        else:
                            print( "Certified unknotted by multiple loops." )
                            return False
                    if foundNewLoop:
                        # Start all over again.
                        break
                    else:
                        print( "Certified unknotted by destroyed loop." )
                        return False
    return
