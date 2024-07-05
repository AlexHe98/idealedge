"""
Tests whether an ideal loop is nontrivially knotted.
"""
from regina import *
from idealedge import decomposeAlong, isSphere
from loop import IdealLoop
try:
    # The multiprocessing package doesn't work with the standard Windows
    # build for Regina.
    from multiprocessing import Process, Pipe
except ModuleNotFoundError:
    warning = "Warning: Proceeding without access to multiprocessing."
    print( "+-" + "-"*len(warning) + "-+" )
    print( "| {} |".format(warning) )
    print( "+-" + "-"*len(warning) + "-+" )
    _serial = True
else:
    from time import sleep
    _serial = False


def knownHyperbolic(loop):
    """
    Is the given ideal loop known to represent a hyperbolic knot?

    The given loop must be embedded in a triangulation of the 3-sphere. Under
    this assumption, if this routine returns True, then the loop is
    guaranteed to be a hyperbolic knot, and is therefore guaranteed to be a
    nontrivial prime knot; otherwise, if this routine returns False, then we
    have no guarantee about whether or not the loop is a hyperbolic knot.
    """
    drilled = loop.drill()
    spt = SnapPeaTriangulation(drilled)
    probablyHyperbolic = False
    attempts = 0
    while True:
        attempts += 1
        if ( spt.solutionType() == SolutionType.geometric_solution or
                spt.solutionType() == SolutionType.nongeometric_solution ):
            probablyHyperbolic = True
            break
        elif attempts < 4:  # Hard-coded limit on the number of attempts.
            # Try again.
            spt.randomise()
        else:
            break
    return ( probablyHyperbolic and spt.hasStrictAngleStructure() )


def isKnotted( loop, tracker=None ):
    """
    Is the given ideal loop nontrivially knotted?

    This routine can be run with an optional DecompositionTracker. The
    intended use case is when a larger decomposition routine needs to track
    progress while running isKnotted() as a subroutine. Thus, this routine
    assumes that tracker.start() has already been called, and it is
    guaranteed that this routine will never call tracker.finish().
    """
    drilled = loop.drill()
    if _serial:
        return _isKnottedSerial( drilled, tracker )
    else:
        return _isKnottedParallel( drilled, tracker )


def _isKnottedParallel( drilled, tracker ):
    # Try enumerating covers on the fundamental group.
    if tracker is not None:
        beforeReport = "Attempting to enumerate covers of index 2 to 6."
        tracker.report(beforeReport)
    isoSig = drilled.isoSig()
    gp = drilled.group()
    for index in range(2,7):
        if tracker is not None:
            tracker.reportIfStalled()

        # If there is more than one cover, then the given triangulation
        # cannot be an ideal solid torus.
        if _tooManyCovers( gp, index, tracker ):
            return True

    # We can still try to enumerate covers of index 7, and if that fails we
    # can resort to solid torus recognition. We use multiprocessing to:
    #   (1) run both computations in parallel; and
    #   (2) facilitate early termination.
    if tracker is not None:
        beforeReport = ( "Enumerating covers of index 7, and " +
                "simultaneously running solid torus recognition." )
        tracker.report(beforeReport)
    coversReceiver, coversSender = Pipe(False)
    notSolidTorusReceiver, notSolidTorusSender = Pipe(False)
    coversProcess = Process(
            target=_runCoversEnumeration, args=( isoSig, 7, coversSender ) )
    notSolidTorusProcess = Process(
            target=_notSolidTorus, args=( isoSig, notSolidTorusSender ) )
    coversProcess.start()
    notSolidTorusProcess.start()
    coversRunning = True
    while True:
        sleep(0.01)
        if tracker is not None:
            try:
                tracker.reportIfStalled()
            except TimeoutError as timeout:
                # Terminate child processes before timing out.
                notSolidTorusProcess.terminate()
                notSolidTorusProcess.join()
                coversProcess.terminate()
                coversProcess.join()
                raise timeout

        # Have we finished enumerating covers of index 7?
        if coversProcess.is_alive():
            # The enumeration is ongoing, but maybe we can terminate early.
            if coversReceiver.poll():
                covers = coversReceiver.recv()
                if covers > 1:
                    coversProcess.terminate()
                    notSolidTorusProcess.terminate()
                    msg = "There are at least {} covers of index 7.".format(
                            covers )
                    isNontrivial = True
                    break
        elif coversRunning:
            coversRunning = False

            # We must have an exact count of the number of covers.
            while coversReceiver.poll():
                covers = coversReceiver.recv()
            if covers != 1:
                notSolidTorusProcess.terminate()
                msg = "Found {} covers of index 7.".format(covers)
                isNontrivial = True
                break

        # Have we finished deciding whether we have the solid torus?
        if not notSolidTorusProcess.is_alive():
            # We must have received a conclusive answer.
            coversProcess.terminate()
            isNontrivial = notSolidTorusReceiver.recv()
            if isNontrivial:
                msg = "Not a solid torus!"
            else:
                msg = "Solid torus!"
            break
    coversProcess.join()
    notSolidTorusProcess.join()
    if tracker is not None:
        tracker.report( None, msg )
    return isNontrivial


def _tooManyCovers( gp, index, tracker ):
    covers = gp.enumerateCovers( index, lambda h: None )
    if tracker is not None and covers != 1:
        afterReport = "Found {} covers of index {}.".format(
                covers, index )
        tracker.report( None, afterReport )
    return ( covers != 1 )


def _runCoversEnumeration( isoSig, index, sender ):
    gp = Triangulation3.fromIsoSig(isoSig).group()
    data = { "covers": 0, "early": False, "sender": sender }
    def handleNewCover(c):
        data["covers"] += 1
        if data["covers"] > 1 and not data["early"]:
            # Notify the parent process that we have already found too many
            # covers, so we can choose to terminate early.
            data["early"] = True
            data["sender"].send( data["covers"] )
        return None
    covers = gp.enumerateCovers( index, handleNewCover )
    sender.send(covers)
    return


def _notSolidTorus( isoSig, sender ):
    drilled = Triangulation3.fromIsoSig(isoSig)
    sender.send( not drilled.isSolidTorus() )
    return


def _isKnottedSerial( drilled, tracker ):
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

        # If there is more than one cover, then the given triangulation
        # cannot be an ideal solid torus.
        if _tooManyCovers( gp, index, tracker ):
            return True

    # Our last resort is to run solid torus recognition directly. Since we
    # survived to this point, the given drilled triangulation is "probably"
    # an ideal solid torus, so this will hopefully terminate quickly.
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


def surgery0(oldLoop):
    """
    Constructs a new ideal loop given by 0/1 Dehn surgery on the given ideal
    loop.

    If the oldLoop is embedded in the 3-sphere, then it is unknotted if and
    only if there is an embedded 2-sphere intersecting the new ideal loop in
    exactly one point.

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
