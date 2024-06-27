"""
Decompose knots into prime knots.
"""
from sys import stdout
from timeit import default_timer
from regina import *
from idealedge import decomposeAlong, isSphere
from loop import IdealLoop
from knotted import isKnotted, knownHyperbolic
try:
    # The multiprocessing package doesn't work with the standard Windows
    # build for Regina.
    from multiprocessing import Process, Pipe
except ModuleNotFoundError:
    _serial = True
else:
    from time import sleep
    _serial = False


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


def reversePinch( knotComplement, packet=None ):
    """
    Builds an ideal loop representing the same knot as the given ideal
    triangulation.

    Warning:
    --> This routine currently uses fast heuristics to attempt to construct
        the desired triangulation, and is not guaranteed to terminate.

    Returns:
        The constructed ideal loop.
    """
    # Triangulate the exterior with boundary edges appearing as the meridian
    # and longitude. The last step is not guaranteed to terminate in theory,
    # but it should be fine in practice.
    knotComplement.minimiseVertices()
    knotComplement.intelligentSimplify()
    knotComplement.intelligentSimplify()
    knotComplement.idealToFinite()
    knotComplement.intelligentSimplify()
    knotComplement.intelligentSimplify()
    mer, lon = knotComplement.meridianLongitude()

    # Get a tetrahedron index and edge number for the longitude, so that we
    # can remember its location after closing up the boundary.
    emb = lon.embedding(0)
    tet = emb.tetrahedron()
    edgeNum = emb.face()

    # Close up the boundary and build the IdealLoop.
    layer = knotComplement.layerOn(mer)
    layer.join( 0, layer, Perm4(0,1) )
    idealEdge = tet.edge(edgeNum)
    loop = IdealLoop( [idealEdge] )
    loop.simplify()
    loop.simplify()
    if packet is not None:
        child = embeddedLoopPacket(loop)
        child.setLabel( packet.adornedLabel(
            "Embedded as edge {}".format( idealEdge.index() ) ) )
        packet.insertChildLast(child)
    return loop


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
    if insertAsChild and isinstance( knot, PacketOfLink ):
        packet = knot
    else:
        packet = None
    knotComplement = knot.complement()
    return reversePinch( knotComplement, packet )


def _runEmbed( knotSig, sender ):
    loop = embedInTriangulation( Link.fromKnotSig(knotSig) )
    sender.send( loop.lightweightDescription() )
    return


def _embedParallel( knot, tracker ):
    receiver, sender = Pipe(False)
    embedProcess = Process(
            target=_runEmbed, args=( knot.knotSig(), sender ) )
    embedProcess.start()
    while True:
        sleep(0.01)
        if tracker is not None:
            try:
                tracker.reportIfStalled()
            except TimeoutError as timeout:
                # Terminate child process before timing out.
                embedProcess.terminate()
                embedProcess.join()
                raise timeout

        # Have we finished the embedding the knot as an ideal loop?
        if not embedProcess.is_alive():
            embedProcess.terminate()
            embedProcess.join()
            loop = IdealLoop()
            loop.setFromLightweight( *receiver.recv() )
            return loop
    return


def _enumerateSerial( oldLoop, tracker ):
    tri = oldLoop.triangulation()
    enumeration = TreeEnumeration( tri, NS_QUAD )
    while True:
        if tracker.hasStalled():
            # We have spent a comparatively long time on the current
            # triangulation, so it might be worthwhile to try harder to
            # simplify this triangulation, and to restart the surface
            # enumeration on a smaller triangulation.
            tracker.report( None, "Try to simplify." )
            drilled = SnapPeaTriangulation( oldLoop.drill() )
            drilled.randomise()
            simpLoop = reversePinch(
                    Triangulation3.fromIsoSig( drilled.isoSig() ) )
            if simpLoop.triangulation().size() < tri.size():
                oldLoop.setFromLoop( simpLoop, False )
                tri = oldLoop.triangulation()
                enumeration = TreeEnumeration( tri, NS_QUAD )
                beforeReport = "Simplified to {} tetrahedra.".format(
                        tri.size() )
                tracker.report(beforeReport)
            else:
                beforeReport = ( "Could not simplify. " +
                        "Continuing with current triangulation." )
                tracker.report(beforeReport)
        tracker.newSearch()

        # Get the next 2-sphere.
        if enumeration.next():
            sphere = enumeration.buildSurface()
            if not isSphere(sphere):
                continue
        else:
            # No suitable 2-sphere means oldLoop is prime.
            # Return ( foundPrime, newLoops ).
            return ( True, None )

        # We only want 2-spheres that intersect the oldLoop in either exactly
        # 0 points or exactly 2 points, since crushing such a 2-sphere has
        # one of the following effects:
        #   --> simplifies the triangulation containing the ideal loop;
        #   --> decomposes the oldLoop into two simpler knots; or
        #   --> (if oldLoop is unknotted) destroys all traces of the loop.
        wt = oldLoop.weight(sphere)
        if wt != 0 and wt != 2:
            continue
        decomposed = decomposeAlong( sphere, [oldLoop] )
        newLoops = []
        for decomposedLoops in decomposed:
            if decomposedLoops:
                # We are guaranteed to have len(decomposedLoops) == 1.
                newLoops.append( decomposedLoops[0] )

        # Return ( foundPrime, newLoops ).
        return ( False, newLoops )


def _perpetualSimplify( isoSig, size, sender ):
    tri = Triangulation3.fromIsoSig(isoSig)
    attempts = 0
    while True:
        attempts += 1
        tri.subdivide()
        tri.intelligentSimplify()
        tri.minimiseVertices()
        tri.intelligentSimplify()
        loop = reversePinch(tri)
        if loop.triangulation().size() < size:
            sender.send( ( loop.lightweightDescription(), attempts ) )
            size = loop.triangulation().size()
            tri = loop.drill()
            attempts = 0
#    spt = SnapPeaTriangulation( Triangulation3.fromIsoSig(isoSig) )
#    attempts = 0
#    while True:
#        attempts += 1
#        spt.randomise()
#        loop = reversePinch( Triangulation3.fromIsoSig( spt.isoSig() ) )
#        if loop.triangulation().size() < size:
#            sender.send( ( loop.lightweightDescription(), attempts ) )
#            size = loop.triangulation().size()
#            spt = SnapPeaTriangulation( loop.drill() )
#            attempts = 0


def _enumerateParallel( oldLoop, tracker ):
    tri = oldLoop.triangulation()
    enumeration = TreeEnumeration( tri, NS_QUAD )
    isoSig = oldLoop.drill().isoSig()
    size = tri.size()
    receiver, sender = Pipe(False)
    simplifyProcess = Process(
            target=_perpetualSimplify, args=( isoSig, size, sender ) )
    simplifyProcess.start()
    while True:
        desc = None
        while receiver.poll():
            desc, attempts = receiver.recv()
        if desc is not None:
            # Successfully simplified! Restart enumeration.
            oldLoop.setFromLightweight( *desc )
            tri = oldLoop.triangulation()
            enumeration = TreeEnumeration( tri, NS_QUAD )
            beforeReport = "Simplified to {} ".format( tri.size() )
            beforeReport += "tetrahedra after "
            if attempts == 1:
                beforeReport += "1 attempt."
            else:
                beforeReport += "{} attempts.".format(attempts)
            try:
                tracker.report(beforeReport)
            except TimeoutError as timeout:
                # Terminate child process before timing out.
                simplifyProcess.terminate()
                simplifyProcess.join()
                raise timeout
        try:
            tracker.newSearch()
        except TimeoutError as timeout:
            # Terminate child process before timing out.
            simplifyProcess.terminate()
            simplifyProcess.join()
            raise timeout

        # Get the next 2-sphere.
        if enumeration.next():
            sphere = enumeration.buildSurface()
            if not isSphere(sphere):
                continue
        else:
            # No suitable 2-sphere means oldLoop is prime.
            # Return ( foundPrime, newLoops ).
            simplifyProcess.terminate()
            simplifyProcess.join()
            return ( True, None )

        # We only want 2-spheres that intersect the oldLoop in either exactly
        # 0 points or exactly 2 points, since crushing such a 2-sphere has
        # one of the following effects:
        #   --> simplifies the triangulation containing the ideal loop;
        #   --> decomposes the oldLoop into two simpler knots; or
        #   --> (if oldLoop is unknotted) destroys all traces of the loop.
        wt = oldLoop.weight(sphere)
        if wt != 0 and wt != 2:
            continue
        decomposed = decomposeAlong( sphere, [oldLoop] )
        newLoops = []
        for decomposedLoops in decomposed:
            if decomposedLoops:
                # We are guaranteed to have len(decomposedLoops) == 1.
                newLoops.append( decomposedLoops[0] )

        # Return ( foundPrime, newLoops ).
        simplifyProcess.terminate()
        simplifyProcess.join()
        return ( False, newLoops )


def decompose( knot, tracker=False, insertAsChild=False ):
    """
    Decomposes the given knot into prime pieces, represented as 3-spheres
    in which the prime knots are embedded as ideal loops.

    The given knot is allowed to be encoded in various ways:
    --> It could be an instance of IdealLoop, in which case it is assumed
        that the triangulation containing this loop is a 3-sphere.
    --> It could be an instance of Regina's Edge3, in which case it is
        assumed that the endpoints of this edge are identified, and that the
        triangulation containing this edge is a 3-sphere.
    --> It could be an instance of Regina's Link or PacketOfLink, in which
        case it is assumed that this link has exactly one component.

    If tracker is an instance of DecompositionTracker, then this routine will
    use this given tracker to track the progress of the decomposition
    computation; if the tracker has the verbose option switched on, then this
    routine will also use the tracker to print regular progress reports.
    Otherwise, the routine will create its own DecompositionTracker, and the
    tracker parameter should be either True or False depending on whether the
    newly-created tracker should have the verbose option switched on.

    If tracker is True or False, then the tracker created by this routine
    will have the timeout feature switched off. Thus, the only way to use the
    timeout feature with this routine is to explicitly supply a tracker.

    If timeout is None (the default), then this routine continues attempting
    to decompose the given knot, regardless of how long this takes.
    Otherwise, timeout should be a positive integer indicating the number of
    seconds after which this routine should attempt to terminate early,
    without completing the decomposition. Such early termination might not
    occur immediately once the allotted time has elapsed, and indeed it might
    not occur at all. Whenever early termination does occur, this routine
    will raise TimeoutError.

    If insertAsChild is True and the given knot is an instance of
    PacketOfLink, then this routine will insert the results of the
    computation as descendents of the given knot packet. This feature is also
    switched off by default.
    """
    if isinstance( tracker, DecompositionTracker ):
        verbose = tracker.isVerbose()
    else:
        verbose = bool(tracker)
        tracker = DecompositionTracker(verbose)
    tracker.start()

    # Build the IdealLoop on which we perform the decomposition computation.
    # Make sure to create clones so as not to directly modify the input.
    if isinstance( knot, IdealLoop ):
        loop = knot.clone()
    elif isinstance( knot, Edge3 ):
        loop = IdealLoop( [knot] ).clone()
    else:
        if verbose:
            beforeReport = "Embedding knot as an ideal loop."
            tracker.report(beforeReport)
        if _serial:
            loop = embedInTriangulation(knot)
        else:
            loop = _embedParallel( knot, tracker )

    # Do the decompositon.
    primes = []
    toProcess = [loop]
    while toProcess:
        # INVARIANT:
        #   At this point, the following are guaranteed to hold:
        #   --> Each element of toProcess is an ideal loop forming a knot.
        #   --> Each element of primes is an ideal loop forming a nontrivial
        #       prime knot.
        #   --> The input knot is given by composing all of the knots
        #       represented in toProcess and primes.
        oldLoop = toProcess.pop()
        tracker.newTri( oldLoop.triangulation().size() )
        if knownHyperbolic(oldLoop):
            # Hyperbolic knots are nontrivial and prime.
            primes.append(oldLoop)
            tracker.foundHyperbolic()
            continue

        # Search for a suitable quadrilateral vertex normal 2-sphere to
        # crush. If no such 2-sphere exists, then the oldLoop is prime.
        # Otherwise, crushing this 2-sphere decomposes the oldLoop into a
        # collection of simpler newLoops.
        if _serial:
            foundPrime, newLoops = _enumerateSerial( oldLoop, tracker )
        else:
            foundPrime, newLoops = _enumerateParallel( oldLoop, tracker )
        if foundPrime:
            # We only care about the case where the prime is nontrivial.
            tracker.unknownPrime()
            isNontrivial = isKnotted( oldLoop, tracker )
            if isNontrivial:
                primes.append(oldLoop)
            tracker.knownPrime(isNontrivial)
        else:
            toProcess.extend(newLoops)
            if verbose:
                tracker.report()

    # Output some auxiliary information before returning the list of primes.
    tracker.finish()
    if verbose:
        msg = tracker.report()
    if insertAsChild and isinstance( knot, PacketOfLink ):
        if verbose:
            container = Text( tracker.log() )
            container.setLabel( "Primes ({})".format(msg) )
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
    A progress tracker for knot decomposition.

    In detail, this tracker provides the following functionality:
    --> Times the tracked knot decomposition computation.
    --> Prints progress reports (either upon request, or upon being notified
        of a significant event).
    --> Tracks whether the computation has stalled, meaning that the number
        of seconds since the most recent event has exceeded some set value.
    --> Provides a timeout option, which raises TimeoutError if one of the
        following occurs after some allotted number of seconds has elapsed:
        (a) This tracker is notified of a significant event.
        (b) This tracker is asked to check whether the computation has
            stalled.

    This tracker recognises the following significant events:
    --> The computation started.
    --> The computation finished.
    --> A progress report was printed.
    --> The computation has begun processing a new triangulation.
    --> The computation has begun a new search for a quadrilateral vertex
        normal surface.
    --> The computation has found a prime knot, but it has not yet
        established whether this prime knot is nontrivially knotted.
    --> The computation has certified whether a prime knot is nontrivially
        knotted.
    """
    def __init__( self, verbose=False, timeout=None, stallInterval=5 ):
        """
        Creates a new DecompositionTracker.

        If verbose is True, then this tracker will automatically print
        progress reports to standard output whenever it is notified of
        significant events; this feature is switched off by default.
        Regardless of whether this feature is switched on or off, it will
        always be possible to manually request a progress report.

        If timeout is None (the default), then the timeout feature will be
        switched off. Otherwise, timeout should be a positive number
        indicating the number of seconds after which the tracked computation
        should be timed out.

        This tracker will consider the tracked knot decomposition computation
        to have stalled if the number of seconds since the last event exceeds
        the given stallInterval.
        """
        self._verbose = verbose
        self._timeout = timeout
        self._indent = "    "
        self._template = "Time: {:.6f}. Searches: {}. Primes: {}. #Tri: {}."
        self._stallInterval = stallInterval
        self._numPrimes = 0
        self._numTri = 0
        self._searches = 0
        self._log = ""
        self._startTime = None
        self._previousEventTime = None
        self._finishTime = None
        return

    def isVerbose(self):
        """
        Is the verbose option switched on for this tracker?
        """
        return self._verbose

    def start(self):
        """
        Starts the timer on the knot decomposition computation that is
        tracked by this tracker.

        This routine must only be called once.
        """
        if self._startTime is not None:
            raise RuntimeError( "Timer already started!" )
        self._startTime = default_timer()
        self._previousEventTime = self._startTime
        return

    def finish(self):
        """
        Informs this tracker that the knot decomposition computation has
        finished.

        This routine must only be called after start() has been called. This
        routine may be called more than once, but calls after the first time
        will do nothing.
        """
        if self._startTime is None:
            raise RuntimeError( "Timer hasn't started yet!" )
        if self._finishTime is not None:
            return
        self._finishTime = default_timer()
        return

    def elapsed(self):
        """
        Returns the total time elapsed during the tracked computation.

        This routine must never be called before start() has been called.
        """
        if self._finishTime is None:
            return default_timer() - self._startTime
        return self._finishTime - self._startTime

    def checkTimeout(self):
        """
        Checks whether the tracked computation should be timed out, and if so
        raises TimeoutError.

        This routine does nothing if the timeout option is switched off or
        the allotted number of seconds has not yet elapsed.
        """
        if self._timeout is not None and self.elapsed() > self._timeout:
            self.finish()
            msg = "Decomposition timed out after {:.6f} seconds.".format(
                    self.elapsed() )
            raise TimeoutError(msg)
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
        self._previousEventTime = time
        msg = self._template.format( time - self._startTime,
                self._searches, self._numPrimes, self._numTri )
        self._printMessage( self._indent + msg )
        return msg

    def report( self, before=None, after=None ):
        """
        Prints and returns a progress report.

        This report may be optionally augmented with messages to appear
        immediately before and/or after the standard progress report.

        This routine raises TimeoutError if it detects that the tracked
        computation should be timed out.

        This routine must never be called before start() has been called.
        """
        if self._finishTime is None:
            time = default_timer()
        else:
            time = self._finishTime
        if before is not None:
            self._printMessage(before)
        rep = self._reportImpl(time)
        if after is not None:
            self._printMessage(after)
        self.checkTimeout()
        return rep

    def _newEvent( self, before=None, after=None ):
        if self._verbose:
            return self.report( before, after )
        self.checkTimeout()
        self._previousEventTime = default_timer()
        return None

    def hasStalled(self):
        """
        Has the tracked computation stalled?

        This routine raises TimeoutError if it detects that the tracked
        computation should be timed out.
        """
        self.checkTimeout()
        if self._finishTime is None:
            return ( default_timer() - self._previousEventTime >
                    self._stallInterval )
        return False

    def _getTimeIfStalled(self):
        self.checkTimeout()
        if self._finishTime is None:
            time = default_timer()
            if time - self._previousEventTime > self._stallInterval:
                return time
        return None

    def reportIfStalled(self):
        """
        Prints and returns a progress report if the tracked computation has
        stalled.

        This routine returns None if the computation is finished, or if the
        computation is still going but has not stalled.

        This routine raises TimeoutError if it detects that the tracked
        computation should be timed out.

        This routine must never be called before start() has been called.
        """
        time = self._getTimeIfStalled()
        if time is not None:
            return self._reportImpl(time)
        return None

    def _newEventIfStalled(self):
        time = self._getTimeIfStalled()
        if time is not None:
            return self._newEvent()
        return None

    def newTri( self, size ):
        """
        Informs this tracker that the tracked computation has started
        processing a new triangulation of the given size.

        If this tracker is verbose, then this routine will automatically
        print a progress report.

        This routine raises TimeoutError if it detects that the tracked
        computation should be timed out.

        This routine must never be called before start() has been called.
        """
        self._numTri += 1
        beforeReport = "Edge-ideal: "
        if size == 1:
            beforeReport += "1 tetrahedron."
        else:
            beforeReport += "{} tetrahedra.".format(size)
        self._newEvent(beforeReport)
        return

    def newSearch(self):
        """
        Informs this tracker that the tracked computation has started a new
        search for a quadrilateral vertex normal surface.

        If this tracker is verbose and the tracked computation has stalled,
        then this routine will automatically print a progress report.

        This routine raises TimeoutError if it detects that the tracked
        computation should be timed out.

        This routine must never be called before start() has been called.
        """
        self._searches += 1
        self._newEventIfStalled()
        return

    def unknownPrime(self):
        """
        Informs this tracker that the tracked computation has found a prime
        knot, but it is not yet known whether this prime is nontrivial.

        If this tracker is verbose, then this routine will automatically
        print a progress report.

        This routine raises TimeoutError if it detects that the tracked
        computation should be timed out.

        This routine must never be called before start() has been called.
        """
        afterReport = "Found a prime knot! Is it nontrivial?"
        self._newEvent( None, afterReport )
        return

    def knownPrime( self, isNontrivial ):
        """
        Informs this tracker that the tracked computation has certified
        whether a prime knot is nontrivially knotted.

        If this tracker is verbose, then this routine will automatically
        print a progress report.

        This routine raises TimeoutError if it detects that the tracked
        computation should be timed out.

        This routine must never be called before start() has been called.
        """
        if isNontrivial:
            self._numPrimes += 1
            beforeReport = "The prime knot is nontrivial!"
        else:
            beforeReport = "The prime knot is the unknot."
        self._newEvent(beforeReport)
        return

    def foundHyperbolic(self):
        """
        Informs this tracker that the tracked computation has found a
        hyperbolic (and hence nontrivial prime) knot.

        If this tracker is verbose, then this routine will automatically
        print a progress report.

        This routine raises TimeoutError if it detects that the tracked
        computation should be timed out.

        This routine must never be called before start() has been called.
        """
        self._numPrimes += 1
        afterReport = "Found a hyperbolic knot!"
        self._newEvent( None, afterReport )
        return
