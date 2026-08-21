"""
Tests whether an edge-ideal triangulation represents a nontrivial knot.
"""
from regina import *
from loop import IdealLoop
from triloops import EdgeIdealTriangulation
from retriangulate.insert import layerOn
from aux.surface import isSphere
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


def knownHyperbolic(edgeIdealTri):
    """
    Is the given edge-ideal triangulation known to represent a hyperbolic
    knot?

    If this routine returns True, then the edge-ideal triangulation is
    guaranteed to represent a hyperbolic knot (and hence a nontrivial prime
    knot); otherwise, if this routine returns False, then we have no
    guarantee about whether or not the knot is hyperbolic.

    Pre-condition:
    --> edgeIdealTri is a 3-sphere containing exactly one ideal loop.
    """
    drilled = edgeIdealTri.drill()
    spt = SnapPeaTriangulation(drilled)
    probablyHyperbolic = False
    attempts = 0
    while True:
        attempts += 1
        sol = spt.solutionType()
        try:
            # Introduced in Regina 7.4:
            geom = SnapPeaTriangulation.Solution.Geometric
            nong = SnapPeaTriangulation.Solution.Nongeometric
        except AttributeError:
            # For backwards compatibility with Regina 7.3 and earlier (but
            # this usage is deprecated as of Regina 7.4):
            geom = SnapPeaTriangulation.geometric_solution
            nong = SnapPeaTriangulation.nongeometric_solution
        if ( sol == geom or sol == nong ):
            probablyHyperbolic = True
            break
        elif attempts < 4:  # Hard-coded limit on the number of attempts.
            # Try again.
            spt.randomise()
        else:
            break
    return ( probablyHyperbolic and spt.hasStrictAngleStructure() )


def isKnotted( edgeIdealTri, tracker=None ):
    """
    Does the given edge-ideal triangulation represent a nontrivial knot?

    This routine can be run with an optional KnotDecompositionTracker. The
    intended use case is when a larger decomposition routine needs to track
    progress while running isKnotted() as a subroutine. Thus, this routine
    assumes that tracker.start() has already been called, and it is
    guaranteed that this routine will never call tracker.finish().

    Pre-condition:
    --> edgeIdealTri is a 3-sphere containing exactly one ideal loop.
    """
    drilled = edgeIdealTri.drill()
    if _serial:
        return _isKnottedSerial( drilled, tracker )
    else:
        return _isKnottedParallel( drilled, tracker )


def _isKnottedParallel( drilled, tracker ):
    # Try enumerating connected k-sheeted covers of the knot exterior.
    if tracker is not None:
        beforeReport = "Attempting to enumerate covers of index 2 to 6."
        tracker.report(beforeReport)
    isoSig = drilled.isoSig()
    gp = drilled.group()
    for index in range(2,7):
        if tracker is not None:
            tracker.reportIfStalled()

        # The fundamental group of the solid torus is Z. Up to conjugacy,
        # this admits exactly one transitive representation into Sym(index),
        # and the abelianisation of the stabiliser subgroup is Z. Thus, if gp
        # does not share these properties, then the given triangulation
        # cannot be an ideal solid torus.
        if _coversDoNotMatch( gp, index, tracker ):
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

        # Have we managed to certify nontriviality using covers of index 7?
        if coversReceiver.poll():
            # The _runCoversEnumeration() routine only sends data if it has
            # successfully certified that the knot is nontrivial.
            notSolidTorusProcess.terminate()
            coversProcess.terminate()
            msg = "Certified nontrivial using covers of index 7."
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


def _coversDoNotMatch( gp, index, tracker ):
    covers = gp.enumerateCovers(index)
    certifiedNontrivial = ( len(covers) != 1 or
            not covers[0].abelianisation().isZ() )
    if tracker is not None and certifiedNontrivial:
        afterReport = ( "Certified nontrivial using " +
                "covers of index {}.".format(index) )
        tracker.report( None, afterReport )
    return certifiedNontrivial


def _runCoversEnumeration( isoSig, index, sender ):
    gp = Triangulation3.fromIsoSig(isoSig).group()
    data = { "sender": sender,
            "covers": 0,
            "certifiedNontrivial": False }
    def handleNewCover(c):
        # No point doing any extra work if we have already certified that the
        # knot is nontrivial.
        if not data["certifiedNontrivial"]:
            data["covers"] += 1
            if data["covers"] > 1 or not c.abelianisation().isZ():
                # Notify the parent process that we have certified that the
                # knot is nontrivial.
                data["certifiedNontrivial"] = True
                data["sender"].send(True)
        return None
    gp.enumerateCovers( index, handleNewCover )
    return


def _notSolidTorus( isoSig, sender ):
    drilled = Triangulation3.fromIsoSig(isoSig)
    sender.send( not drilled.isSolidTorus() )
    return


def _isKnottedSerial( drilled, tracker ):
    # Try enumerating connected k-sheeted covers of the knot exterior.
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

        # The fundamental group of the solid torus is Z. Up to conjugacy,
        # this admits exactly one transitive representation into Sym(index),
        # and the abelianisation of the stabiliser subgroup is Z. Thus, if gp
        # does not share these properties, then the given triangulation
        # cannot be an ideal solid torus.
        if _coversDoNotMatch( gp, index, tracker ):
            return True

    # Our last resort is to run solid torus recognition directly. Since we
    # survived to this point, the given drilled triangulation is "probably"
    # an ideal solid torus, so this will hopefully terminate quickly.
    drilled.idealToFinite()
    #NOTE Triangulation3.simplify() was introduced in Regina 7.4. In older
    #       versions of Regina, equivalent functionality was provided by
    #       Triangulation3.intelligentSimplify().
    drilled.simplify()
    drilled.simplify()
    if tracker is not None:
        beforeReport = "Resorting to solid torus recognition.\n"
        beforeReport += "Truncated: {} tetrahedra.".format( drilled.size() )
        tracker.report(beforeReport)
    isNontrivial = not drilled.isSolidTorus()
    if tracker is not None:
        tracker.report()
    return isNontrivial


#TODO Do we still need this?
def surgery0(oldEdgeIdealTri):
    """
    Constructs a new edge-ideal triangulation given by 0/1 Dehn surgery on
    the given edge-ideal triangulation.

    The given oldEdgeIdealTri must contain exactly one ideal loop. If the
    ambient triangulation is a 3-sphere, then the loop is unknotted if and
    only if the returned triangulation contains an embedded 2-sphere
    intersecting the newly-constructed ideal loop in exactly one point.

    This routine might raise BoundsDisc.

    Warning:
    --> This routine currently uses fast heuristics to attempt to construct
        the desired triangulation, and is not guaranteed to terminate.

    Returns:
        The newly-constructed edge-ideal triangulation.
    """
    # Triangulate the exterior with boundary edges appearing as the meridian
    # and longitude. The last step is not guaranteed to terminate in theory,
    # but it should be fine in practice.
    tri = oldEdgeIdealTri.drill()
    tri.idealToFinite()
    #NOTE Triangulation3.simplify() was introduced in Regina 7.4. In older
    #       versions of Regina, equivalent functionality was provided by
    #       Triangulation3.intelligentSimplify().
    tri.simplify()
    tri.simplify()
    mer, lon = tri.meridianLongitude()

    # Get a tetrahedron index and edge number for the meridian, so that we
    # can remember its location after closing up the boundary.
    emb = mer.embedding(0)
    tet = emb.tetrahedron()
    edgeNum = emb.face()

    # Close up the boundary and build the new EdgeIdealTriangulation.
    #
    # newEdgeIdealTri.simplify() might raise BoundsDisc.
    layer = layerOn(lon)
    layer.join( 0, layer, Perm4(0,1) )
    idealEdge = tet.edge(edgeNum)
    newEdgeIdealTri = EdgeIdealTriangulation( [ IdealLoop( [idealEdge] ) ] )
    newEdgeIdealTri.simplify()
    newEdgeIdealTri.simplify()
    return newEdgeIdealTri
