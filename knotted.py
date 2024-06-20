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

    This routine can be run with an optional DecompositionTracker. The
    intended use case is when a larger decomposition routine needs to track
    progress while running isKnotted() as a subroutine. Thus, this routine
    assumes that tracker.start() has already been called, and it is
    guaranteed that this routine will never call tracker.finish().
    """
    # If we can certify hyperbolicity, then the loop must be knotted.
    drilled = loop.drill()
    if tracker is not None:
        msg = "Drilled: {} tetrahedra.\n".format( drilled.size() )
        msg += "Attempting to certify hyperbolic."
        tracker.report(msg)
    spt = SnapPeaTriangulation(drilled)
    probablyHyperbolic = False
    attempts = 0
    while True:
        if tracker is not None:
            tracker.reportIfStalled()
        attempts += 1
        if ( spt.solutionType() == SolutionType.geometric_solution or
                spt.solutionType() == SolutionType.nongeometric_solution ):
            probablyHyperbolic = True
            break
        elif attempts < 4:
            # Try again.
            spt.randomise()
        else:
            break
    if probablyHyperbolic and spt.hasStrictAngleStructure():
        # Certified hyperbolic.
        if tracker is not None:
            tracker.report()
        return True

    #
    #TODO
    core = surgery0(loop)
    while True:
        # INVARIANT:
        #   It is guaranteed that the original loop is unknotted if and only
        #   if there is a normal 2-sphere intersecting the current core loop
        #   in *exactly* one point.
        tri = core.triangulation()
        if tracker is not None:
            tracker.newTri( tri.size() )

        # Search for a quadrilateral vertex normal 2-sphere that intersects
        # the current core loop in *at most* one point.
        # - If no such 2-sphere exists, then the original loop is
        #   nontrivially knotted.
        # - If we find a 2-sphere intersecting the current core loop in
        #   *exactly* one point, then the original loop is unknotted.
        # - If we find a 2-sphere disjoint from the current core loop, then
        #   we can crush to obtain a new core loop that still satisfies the
        #   invariant.
        enumeration = TreeEnumeration( tri, NS_QUAD )
        while True:
            if tracker is not None:
                tracker.newSearch()

            # Get the next 2-sphere.
            if enumeration.next():
                sphere = enumeration.buildSurface()
                if not isSphere(sphere):
                    continue
            else:
                # As above, no suitable 2-sphere means the original loop is
                # nontrivially knotted.
                return True

            # As above, we only want at most one point of intersection with
            # the core loop.
            wt = core.weight(sphere)
            if wt == 1:
                # The original loop is unknotted.
                return False
            elif wt == 0:
                # Crushing is guaranteed to give us exactly one new component
                # containing an ideal loop, and this new ideal loop will be
                # our new core loop.
                decomposed = decomposeAlong( sphere, [core] )
                for newLoops in decomposed:
                    if newLoops:
                        # We are guaranteed to have len(newLoops) == 1.
                        core = newLoops[0]
                        break
                break
    return
