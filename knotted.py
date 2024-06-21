"""
Uses an ideal loop given by 0/1 Dehn surgery to test whether a knot is
nontrivially knotted.
"""
from regina import *
from idealedge import decomposeAlong, isSphere
from loop import IdealLoop
try:
    # The multiprocessing package doesn't work with the standard Windows
    # build for Regina.
    from multiprocessing import Process, Pipe
except ModuleNotFoundError:
    print( "Warning: Proceeding without access to multiprocessing." )
    _inSeries = True
else:
    _inSeries = False


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
    newLoop.simplify()
    return newLoop


def _countCovers( gp, index ):
    return gp.enumerateCovers( index, print )


def _notSolidTorus():
    #TODO
    return


def _isKnottedInSeries( drilled, tracker ):
    # Try enumerating covers on the fundamental group.
    if tracker is not None:
        beforeReport = "Attempting to enumerate covers of index 2 to 6."
        tracker.report(beforeReport)
    gp = drilled.group()
    for index in range(2,8):
        if tracker is not None:
            if index == 7:
                beforeReport = "Attempting to enumerate covers of index 7."
                tracker.report(beforeReport)
            else:
                tracker.reportIfStalled()
        covers = _countCovers( gp, index )
        if covers != 1:
            # The unknot has only one cover, so the given loop must be
            # nontrivially knotted.
            if tracker is not None:
                afterReport = "Found {} covers of index {}!".format(
                        covers, index )
                tracker.report( None, afterReport )
            return True

    # If we survive to this point, then the given drilled triangulation is
    # probably not the complement of a hyperbolic knot, and algebra suggests
    # that it is likely to unknotted. Anyway, our last resort is to run solid
    # torus recognition directly.
    drilled.idealToFinite()
    drilled.intelligentSimplify()
    drilled.intelligentSimplify()
    if tracker is not None:
        beforeReport = "Resorting to solid torus recognition.\n"
        beforeReport += "Truncated: {} tetrahedra.".format( drilled.size() )
        tracker.report(beforeReport)
    isNontrivial = not drilled.isSolidTorus()
    if tracker is not None:
        tracker.report()
    return isNontrivial


def _isKnottedParallel():
    #TODO Try using multiprocessing to enable both parallel computation and
    #   early termination.
    raise NotImplementedError()


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
        beforeReport = "Drilled: {} tetrahedra.\n".format( drilled.size() )
        beforeReport += "Attempting to certify hyperbolic."
        tracker.report(beforeReport)
    spt = SnapPeaTriangulation(drilled)
    probablyHyperbolic = False
    attempts = 0
    while True:
        if tracker is not None:
            # Do this just in case, even though we do not expect these
            # hyperbolic tests to take a long time.
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
            afterReport = "Certified hyperbolic!"
            tracker.report( None, afterReport )
        return True

    # Now try enumerating covers on the fundamental group and/or simply
    # resorting to running solid torus recognition directly.
    if _inSeries:
        return _isKnottedInSeries( drilled, tracker )
    else:
        return _isKnottedParallel( drilled, tracker )
