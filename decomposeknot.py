"""
Decompose knots into prime knots.
"""
from sys import stdout
from timeit import default_timer
from regina import *
from idealedge import decomposeAlong, isSphere
from loop import IdealLoop


def embeddedLoopPacket(loop):
    """
    Returns a packet of the triangulation containing the given loop, with an
    ideal triangulation of the drilled 3-manifold as a child.
    """
    drilled = PacketOfTriangulation3( loop.drill() )
    drilled.setLabel( "Drilled: {}".format( drilled.isoSig() ) )
    packet = PacketOfTriangulation3( loop.triangulation() )
    packet.insertChildLast(drilled)
    return packet


def embedInTriangulation( knot, insertAsChild=False ):
    """
    Embeds the given knot as an ideal loop in a triangulation of the
    3-sphere.

    Warning:
    --> This routine currently uses fast heuristics to attempt to construct
        the desired triangulation, and is not guaranteed to terminate.

    Returns:
        The constructed ideal loop.
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
    tri.intelligentSimplify()
    mer, lon = tri.meridianLongitude()

    # Get a tetrahedron index and edge number for the longitude, so that we
    # can remember its location after closing up the boundary.
    emb = lon.embedding(0)
    tet = emb.tetrahedron()
    edgeNum = emb.face()

    # Close up the boundary and build the IdealLoop.
    layer = tri.layerOn(mer)
    layer.join( 0, layer, Perm4(0,1) )
    idealEdge = tet.edge(edgeNum)
    loop = IdealLoop( [idealEdge] )
    loop.simplify()
    if insertAsChild and isinstance( knot, PacketOfLink ):
        packet = embeddedLoopPacket(loop)
        packet.setLabel( knot.adornedLabel(
            "Embedded as edge {}".format( idealEdge.index() ) ) )
        knot.insertChildLast(packet)
    return loop


def decompose( knot, verbose=False, insertAsChild=False ):
    """
    Decomposes the given knot into prime pieces, represented as 3-spheres
    in which the prime knots are embedded as ideal loops.

    If verbose is True, then this routine will print regular progress
    reports. If insertAsChild is True and the given knot is an instance of
    PacketOfLink, then this routine will insert the results of the
    computation as descendents of the given knot packet. Both of these
    features are switched off by default.
    """
    if verbose:
        tracker = DecompositionTracker()
        tracker.start()
    primes = []
    toProcess = [ embedInTriangulation(knot) ]
    while toProcess:
        # INVARIANT:
        #   At this point, the following are guaranteed to hold:
        #   --> Each element of toProcess is an ideal loop forming a knot.
        #   --> Each element of primes is an ideal loop forming a nontrivial
        #       prime knot.
        #   --> The input knot is given by composing all of the knots
        #       represented in toProcess and primes.
        oldLoop = toProcess.pop()
        tri = oldLoop.triangulation()
        if verbose:
            tracker.newTri( tri.size() )

        # Search for a suitable quadrilateral vertex normal 2-sphere to
        # crush. If no such 2-sphere exists, then the oldLoop is prime.
        enumeration = TreeEnumeration( tri, NS_QUAD )
        while True:
            if verbose:
                tracker.newSearch()

            # Get the next 2-sphere.
            if enumeration.next():
                sphere = enumeration.buildSurface()
                if not isSphere(sphere):
                    continue
            else:
                # No suitable 2-sphere means oldLoop is prime. But we only
                # care about the case where this prime is nontrivial.
                if verbose:
                    drilled = oldLoop.drill()
                    drilled.idealToFinite()
                    drilled.intelligentSimplify()
                    drilled.intelligentSimplify()
                    tracker.unknownPrime( drilled.size() )
                    #TODO Use ideal-edge machinery here.
                    isKnotted = not drilled.isSolidTorus()
                    if isKnotted:
                        primes.append(oldLoop)
                    tracker.knownPrime(isKnotted)
                elif not oldLoop.drill().isSolidTorus():
                    primes.append(oldLoop)
                break

            # We only want 2-spheres that intersect the oldLoop in either
            # exactly 0 points or exactly 2 points, since crushing such a
            # 2-sphere either:
            # - simplifies the triangulation containing the ideal loop;
            # - decomposes the oldLoop into two simpler knots; or
            # - (if oldLoop is unknotted) destroys all traces of the loop.
            wt = oldLoop.weight(sphere)
            if wt != 0 and wt != 2:
                continue
            decomposed = decomposeAlong( sphere, [oldLoop] )
            knots = []
            for newLoops in decomposed:
                if newLoops:
                    # We are guaranteed to have len(newLoops) == 1.
                    knots.append( newLoops[0] )
            for newLoop in knots:
                toProcess.append(newLoop)
            if verbose:
                tracker.report()
            break

    # Output some auxiliary information before returning the list of primes.
    if verbose:
        tracker.report()
    if insertAsChild and isinstance( knot, PacketOfLink ):
        if verbose:
            container = Text( tracker.log() )
        else:
            container = Container()
        container.setLabel("Primes")
        knot.insertChildLast(container)
        for i, primeLoop in enumerate(primes):
            packet = embeddedLoopPacket(primeLoop)
            loopEdgeIndices = list(primeLoop)
            if len(primeLoop) == 1:
                adorn = "Embedded as edge {}".format( loopEdgeIndices[0] )
            else:
                indices = ""
                for ei in loopEdgeIndices[:-1]:
                    indices += ", {}".format(ei)
                adorn = "Embedded as edges {} and {}".format(
                        indices[2:], loopEdgeIndices[-1] )
            packet.setLabel( "Prime knot #{} ({})".format( i, adorn ) )
            container.insertChildLast(packet)
    return primes


class DecompositionTracker:
    """
    A progress tracker for knot decomposition, whose main purpose is to print
    progress reports when running decompose() with the verbose option.
    """
    def __init__( self, stallInterval=5 ):
        """
        Create a new DecompositionTracker with the given stallInterval.

        This DecompositionTracker will consider the tracked knot
        decomposition computation to have stalled if the number of seconds
        since the previous progress report exceeds the stallInterval.
        """
        self._template = ( "    " +
                "Time: {:.6f}. Searches: {}. Primes: {}. #Tri: {}." )
        self._stallInterval = stallInterval
        self._numPrimes = 0
        self._numTri = 0
        self._searches = 0
        self._log = ""
        self._started = False
        return

    def start(self):
        """
        Starts the timer on the knot decomposition computation that is
        tracked by this tracker.

        This routine must only be called once.
        """
        if self._started:
            raise RuntimeError( "Timer already started!" )
        self._started = True
        self._start = default_timer()
        self._prev = self._start
        return

    def log(self):
        """
        Returns a log of all progress reports that have appeared so far.

        The log will be a string that could consist of many lines of text.
        """
        return self._log

    def _printMessage( self, msg ):
        self._log += msg + "\n"
        print(msg)
        stdout.flush()
        return

    def _reportImpl( self, time ):
        self._prev = time
        msg = self._template.format( time - self._start,
                self._searches, self._numPrimes, self._numTri )
        self._printMessage(msg)
        return

    def report(self):
        """
        Prints a progress report.
        """
        self._reportImpl( default_timer() )
        return

    def reportIfStalled(self):
        """
        Prints a progress report if the tracked computation has stalled.
        """
        time = default_timer()
        if time - self._prev > self._stallInterval:
            self._reportImpl(time)
        return

    def newTri( self, size ):
        """
        Informs this tracker that the tracked computation has started
        processing a new triangulation of the given size.

        This routine will also automatically print a progress report.
        """
        self._numTri += 1
        msg = "Edge-ideal: "
        if size == 1:
            msg += "1 tetrahedron."
        else:
            msg += "{} tetrahedra.".format(size)
        self._printMessage(msg)
        self.report()
        return

    def newSearch(self):
        """
        Informs this tracker that the tracked computation has started a new
        search for a quadrilateral vertex normal surface.

        This routine will also automatically print a progress report if the
        tracked computation has stalled.
        """
        self._searches += 1
        self.reportIfStalled()
        return

    def unknownPrime( self, size ):
        """
        Informs this tracker that the tracked computation has found a prime
        knot, but it is not yet known whether this prime is nontrivial.

        The given size should be the number of tetrahedra in an ideal
        triangulation of the found prime knot.

        This routine will also automatically print a progress report.
        """
        self.report()
        msg = "Found prime knot! Is it nontrivial?\nDrilled: "
        if size == 1:
            msg += "1 tetrahedron."
        else:
            msg += "{} tetrahedra.".format(size)
        self._printMessage(msg)
        return

    def knownPrime( self, isKnotted ):
        """
        Informs this tracker that the tracked computation has certified
        whether a prime knot is nontrivially knotted.

        This routine will also automatically print a progress report.
        """
        if isKnotted:
            self._numPrimes += 1
            msg = "Yes, found a nontrivial prime knot!"
        else:
            msg = "No, it's the unknot."
        self.report()
        self._printMessage(msg)
        return
