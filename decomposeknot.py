"""
Decompose knots into prime knots.
"""
from sys import stdout
from timeit import default_timer
from regina import *
from idealedge import decomposeAlong
from loop import IdealLoop, BoundsDisc
from triloops import EdgeIdealTriangulation
from knotted import isKnotted, knownHyperbolic
from embed import loopPacket, reversePinch, embedByFilling, embedFromDiagram
from aux.surface import isSphere
try:
    # The multiprocessing package doesn't work with the standard Windows
    # build for Regina.
    from multiprocessing import Process, Pipe
except ModuleNotFoundError:
    _serial = True
else:
    from time import sleep
    _serial = False


def decompose( knot, tracker=False, insertAsChild=False ):
    """
    Decomposes the given knot into prime summands, represented as edge-ideal
    triangulations.

    In detail, this routine returns a list of EdgeIdealTriangulation objects,
    each of which is a 3-sphere with one of the prime summands embedded as an
    IdealLoop.

    The given knot is allowed to be encoded in various ways:
    --> It could be an instance of EdgeIdealTriangulation, in which case it
        is assumed that this consists of a 3-sphere triangulation containing
        exactly one IdealLoop.
    --> It could be an instance of Regina's Edge3, in which case it is
        assumed that the endpoints of this edge are identified, and that the
        triangulation containing this edge is a 3-sphere.
    --> It could be an instance of Regina's Link or PacketOfLink, in which
        case it is assumed that this link has exactly one component.

    If tracker is an instance of KnotDecompositionTracker, then this routine
    will use this given tracker to track the progress of the decomposition
    computation; if the tracker has the verbose option switched on, then this
    routine will also use the tracker to print regular progress reports.
    Otherwise, the routine will create its own KnotDecompositionTracker, and
    the tracker parameter should be either True or False depending on whether
    the newly-created tracker should have the verbose option switched on.

    If tracker is True or False, then the tracker created by this routine
    will have the timeout feature switched off. Thus, the only way to use the
    timeout feature with this routine is to explicitly supply a tracker.

    An explicitly supplied tracker may or may not already be started. If it
    is already started, then it will be assumed that the tracker is tracking
    a larger computation; hence, this routine will not call the tracker's
    finish() routine. Otherwise, if the tracker is not already started, then
    it will be assumed that the tracker is intended to track only the
    progress of this routine; hence, this routine will call the tracker's
    start() routine before performing the bulk of the computation, and it
    will also call the tracker's finish() routine after it has performed the
    bulk of the computation.

    If insertAsChild is True and the given knot is an instance of
    PacketOfLink, then this routine will insert the results of the
    computation as descendents of the given knot packet. This feature is also
    switched off by default.
    """
    if isinstance( tracker, KnotDecompositionTracker ):
        verbose = tracker.isVerbose()
    else:
        verbose = bool(tracker)
        tracker = KnotDecompositionTracker(verbose)
    if tracker.isStarted():
        needToFinish = False
    else:
        needToFinish = True
        tracker.start()

    # Build the EdgeIdealTriangulation on which we perform the decomposition
    # computation. Make sure to create clones so as not to directly modify
    # the input.
    if isinstance( knot, EdgeIdealTriangulation ):
        edgeIdealTri = knot.clone()
    elif isinstance( knot, Edge3 ):
        edgeIdealTri = EdgeIdealTriangulation( IdealLoop( [knot] ).clone() )
    else:
        if verbose:
            beforeReport = "Knot sig: {}.\n".format( knot.knotSig() )
            beforeReport += "Embedding knot as an ideal loop."
            tracker.report(beforeReport)
        if _serial:
            # In practice, embedFromDiagram(knot) is usually slower than
            # embedByFilling(knot). However, embedFromDiagram(knot) is
            # guaranteed to terminate, so it is the better option if we are
            # not able to use multiprocessing.
            try:
                edgeIdealTri = embedFromDiagram(knot)
            except BoundsDisc:
                # The given knot is unknotted.
                edgeIdealTri = None
        else:
            try:
                edgeIdealTri = _embedParallel( knot, tracker )
            except BoundsDisc:
                # The given knot is unknotted.
                edgeIdealTri = None

    # Do the decompositon.
    primes = []
    if edgeIdealTri is None:
        # The given knot is unknotted.
        if verbose:
            afterReport = "The knot bounds a disc!"
            tracker.report( None, afterReport )
        toProcess = []
    else:
        toProcess = [edgeIdealTri]
    while toProcess:
        # INVARIANT:
        #   At this point, the following are guaranteed to hold:
        #   --> Each element of toProcess is an edge-ideal triangulation
        #       representing a (possibly trivial, possibly composite) knot.
        #   --> Each element of primes is an edge-ideal triangulation
        #       representing a nontrivial prime knot.
        #   --> The input knot is given by composing all of the knots
        #       represented in toProcess and primes.
        oldEdgeIdealTri = toProcess.pop()
        tracker.newEdgeIdealTri(oldEdgeIdealTri)
        if knownHyperbolic(oldEdgeIdealTri):
            # Hyperbolic knots are nontrivial and prime.
            primes.append(oldEdgeIdealTri)
            tracker.foundHyperbolic()
            continue

        # Search for a suitable quadrilateral vertex normal 2-sphere to
        # crush. If no such 2-sphere exists, then the oldEdgeIdealTri is
        # prime. Otherwise, crushing this 2-sphere decomposes the
        # oldEdgeIdealTri into a collection of simpler newEdgeIdealKnots.
        if _serial:
            try:
                newEdgeIdealKnots = _enumerateSerial(
                        oldEdgeIdealTri, tracker )
            except BoundsDisc:
                # The oldEdgeIdealTri is unknotted.
                tracker.knownPrime(False)
                continue
            msg = None
        else:
            try:
                newEdgeIdealKnots, msg = _enumerateParallel(
                        oldEdgeIdealTri, tracker )
            except BoundsDisc:
                # The oldEdgeIdealTri is unknotted.
                tracker.knownPrime(False)
                continue
        if newEdgeIdealKnots is None:
            # The oldEdgeIdealTri is prime! However, we only care about the case
            # where this prime is nontrivial.
            tracker.unknownPrime(msg)
            isNontrivial = isKnotted( oldEdgeIdealTri, tracker )
            if isNontrivial:
                primes.append(oldEdgeIdealTri)
            tracker.knownPrime(isNontrivial)
        else:
            toProcess.extend(newEdgeIdealKnots)
            if verbose:
                tracker.report( None, msg )

    # Output some auxiliary information before returning the list of primes.
    if needToFinish:
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
        for i, primeEdgeIdealTri in enumerate(primes):
            packet = loopPacket(primeEdgeIdealTri)
            loopEdgeIndices = list( primeEdgeIdealTri[0] )
            if len( primeEdgeIdealTri[0] ) == 1:
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


def _embedParallel( knot, tracker ):
    knotSig = knot.knotSig()

    # Run embedByFilling() in a child process.
    fillingReceiver, fillingSender = Pipe(False)
    fillingProcess = Process(
            target=_runFilling, args=( knotSig, fillingSender ) )
    fillingProcess.start()

    # Run embedFromDiagram() in a child process.
    diagramReceiver, diagramSender = Pipe(False)
    diagramProcess = Process(
            target=_runDiagram, args=( knotSig, diagramSender ) )
    diagramProcess.start()

    # The Hare and the Tortoise
    # -------------------------
    # The advantage of embedByFilling() is that it is faster in most cases,
    # whereas the advantage of embedFromDiagram() is that it is guaranteed to
    # terminate. We don't care who wins the race; both will give us a
    # suitable edge-ideal triangulation.
    while True:
        sleep(0.01)
        if tracker is not None:
            try:
                tracker.reportIfStalled()
            except TimeoutError as timeout:
                # Terminate child processes before timing out.
                fillingProcess.terminate()
                diagramProcess.terminate()
                fillingProcess.join()
                diagramProcess.join()
                raise timeout

        # Have we finished building an edge-ideal triangulation for the knot?
        if not fillingProcess.is_alive():
            diagramProcess.terminate()
            fillingProcess.join()
            diagramProcess.join()
            if fillingReceiver.poll():
                edgeIdealTri = EdgeIdealTriangulation.fromBlueprint(
                        *fillingReceiver.recv() )
            else:
                # If fillingProcess terminated without sending information,
                # then the given knot must be unknotted.
                raise BoundsDisc()
            if tracker is not None:
                afterReport = "Built triangulation using 1/0 Dehn surgery."
                tracker.report( None, afterReport )
            return edgeIdealTri
        if not diagramProcess.is_alive():
            fillingProcess.terminate()
            diagramProcess.join()
            fillingProcess.join()
            if diagramReceiver.poll():
                edgeIdealTri = EdgeIdealTriangulation.fromBlueprint(
                        *diagramReceiver.recv() )
            else:
                # If diagramProcess terminated without sending information,
                # then the given knot must be unknotted.
                raise BoundsDisc()
            if tracker is not None:
                afterReport = "Built triangulation from planar diagram."
                tracker.report( None, afterReport )
            return edgeIdealTri
    return


def _runFilling( knotSig, sender ):
    RandomEngine.reseedWithHardware()
    try:
        edgeIdealTri = embedByFilling( Link.fromKnotSig(knotSig) )
    except BoundsDisc:
        # Send nothing if the given knot is unknotted.
        return
    sender.send( edgeIdealTri.blueprint() )
    return


def _runDiagram( knotSig, sender ):
    RandomEngine.reseedWithHardware()
    try:
        edgeIdealTri = embedFromDiagram( Link.fromKnotSig(knotSig) )
    except BoundsDisc:
        # Send nothing if the given knot is unknotted.
        return
    sender.send( edgeIdealTri.blueprint() )
    return


def _enumerateParallel( oldEdgeIdealTri, tracker ):
    # Searching for quadrilateral vertex normal 2-spheres can be very slow.
    # However, if the oldEdgeIdealTri is a composite knot, then in practice
    # we find that we can "usually" find the desired 2-sphere very quickly.
    # Thus, when the enumeration takes a long time for the given
    # oldEdgeIdealTri, it is often helpful to randomise the edge-ideal
    # triangulation and attempt the enumeration on the new triangulation.
    blueprint = oldEdgeIdealTri.blueprint()
    tri = oldEdgeIdealTri.triangulation()

    # Set up a child process to repeatedly randomise the given edge-ideal
    # triangulation, and send the randomised triangulations to another child
    # process that runs alternate enumerations.
    randomiseReceiver, randomiseSender = Pipe(False)
    randomiseProcess = Process( target=_perpetualRandomise,
            args=( blueprint, tri.size(), randomiseSender ) )
    randomiseProcess.start()

    # Set up a child process to run the alternate enumerations.
    alternateReceiver, alternateSender = Pipe(False)
    alternateProcess = Process( target=_indefiniteEnumerate,
            args=( randomiseReceiver, alternateSender ) )
    alternateProcess.start()

    # Run the main enumeration.
    enumeration = TreeEnumeration( tri, NS_QUAD )
    msg = "Main enumeration succeeded."
    while True:
        # Has the randomiseProcess determined that the oldEdgeIdealTri
        # represents an unknot?
        if not randomiseProcess.is_alive():
            # Make sure to clean up child processes before raising BoundsDisc
            # to indicate that the oldEdgeIdealTri is unknotted.
            alternateProcess.terminate()
            randomiseProcess.join()
            alternateProcess.join()
            raise BoundsDisc()

        # Has the alternateProcess given an answer?
        if alternateReceiver.poll():
            # Make sure to clean up child processes before returning the
            # answer from the alternateProcess.
            randomiseProcess.terminate()
            alternateProcess.join()
            randomiseProcess.join()
            newBlueprints, attempts, searches, size =\
                    alternateReceiver.recv()
            msg = "Alternate enumeration succeeded on "
            msg += "{}-tetrahedron triangulation.\n".format(size)
            msg += "(Randomisation attempts: {}. Searches: {}.)".format(
                    attempts, searches )
            if newBlueprints is None:
                # Found a prime!
                return ( None, msg )
            else:
                # Build new edge-ideal triangulations and return them.
                newEdgeIdealKnots = []
                for blueprint in newBlueprints:
                    edgeIdealTri = EdgeIdealTriangulation.fromBlueprint(
                            *blueprint )
                    newEdgeIdealKnots.append(edgeIdealTri)
                return ( newEdgeIdealKnots, msg )

        # Continue with main enumeration (if not timed out).
        try:
            tracker.newSearch()
        except TimeoutError as timeout:
            # Terminate child processes before timing out.
            alternateProcess.terminate()
            randomiseProcess.terminate()
            alternateProcess.join()
            randomiseProcess.join()
            raise timeout

        # Get the next 2-sphere.
        if enumeration.next():
            sphere = enumeration.buildSurface()
            if not isSphere(sphere):
                continue
        else:
            # No suitable 2-sphere means oldEdgeIdealTri is prime.
            # Clean up child processes before returning.
            alternateProcess.terminate()
            randomiseProcess.terminate()
            alternateProcess.join()
            randomiseProcess.join()
            return ( None, msg )

        # We only want 2-spheres that intersect the oldEdgeIdealTri in either
        # exactly 0 points or exactly 2 points, since crushing such a
        # 2-sphere has one of the following effects:
        #   --> simplifies the triangulation containing the ideal loop;
        #   --> decomposes the oldEdgeIdealTri into two simpler knots; or
        #   --> (if oldEdgeIdealTri is unknotted) destroys all traces of the
        #       loop.
        wt = oldEdgeIdealTri.weight(sphere)
        if wt != 0 and wt != 2:
            continue
        # Because we are working with edge-ideal triangulations of knots in
        # the 3-sphere, orbital compression discs and deleted components
        # correspond to either 3-spheres or edge-ideal triangulations of
        # unknots, so we can ignore all of the extra bookkeeping.
        decomposed, numOrbCuts, delComps, inconsistent = decomposeAlong(
                sphere, EdgeIdealTriangulation( [oldEdgeIdealTri] ) )
        newEdgeIdealKnots = []
        for edgeIdealTri in decomposed:
            if isinstance( edgeIdealTri, EdgeIdealTriangulation ):
                try:
                    edgeIdealTri.simplify()
                except BoundsDisc:
                    # We can just throw away this unknot.
                    continue
                else:
                    # We are guaranteed to have len(edgeIdealTri) == 1.
                    newEdgeIdealKnots.append(edgeIdealTri)

        # Clean up child processes before returning.
        alternateProcess.terminate()
        randomiseProcess.terminate()
        alternateProcess.join()
        randomiseProcess.join()
        return ( newEdgeIdealKnots, msg )
    return


def _perpetualRandomise( blueprint, size, sender ):
    RandomEngine.reseedWithHardware()
    edgeIdealTri = EdgeIdealTriangulation.fromBlueprint( *blueprint )
    attempts = 0
    while True:
        attempts += 1
        try:
            edgeIdealTri.randomise()    # Might raise BoundsDisc.
        except BoundsDisc:
            # Use early termination to indicate that the edgeIdealTri
            # represents an unknot.
            return
        if edgeIdealTri.triangulation().size() <= size:
            # Send randomised edgeIdealTri.
            sender.send( ( edgeIdealTri.blueprint(), attempts ) )
    return


def _indefiniteEnumerate( receiver, sender ):
    searches = 0
    while not receiver.poll():
        sleep(0.01)
    blueprint, attempts = receiver.recv()
    edgeIdealTri = EdgeIdealTriangulation.fromBlueprint( *blueprint )
    tri = edgeIdealTri.triangulation()
    enumeration = TreeEnumeration( tri, NS_QUAD )
    while True:
        if searches > 20 and receiver.poll():
            # Restart the enumeration with a new edge-ideal triangulation.
            searches = 0
            blueprint, attempts = receiver.recv()
            edgeIdealTri = EdgeIdealTriangulation.fromBlueprint( *blueprint )
            tri = edgeIdealTri.triangulation()
            enumeration = TreeEnumeration( tri, NS_QUAD )

        # Get the next 2-sphere.
        searches += 1
        if enumeration.next():
            sphere = enumeration.buildSurface()
            if not isSphere(sphere):
                continue
        else:
            # No suitable 2-sphere means the edge-ideal triangulation
            # represents a prime knot.
            sender.send( ( None, attempts, searches, tri.size() ) )
            return

        # We only want 2-spheres that intersect the loop in either exactly 0
        # points or exactly 2 points, since crushing such a 2-sphere has one
        # of the following effects:
        #   --> simplifies the triangulation containing the ideal loop;
        #   --> decomposes the loop into two simpler knots; or
        #   --> (if the loop is unknotted) destroys all traces of the loop.
        wt = edgeIdealTri.weight(sphere)
        if wt != 0 and wt != 2:
            continue
        # Because we are working with edge-ideal triangulations of knots in
        # the 3-sphere, orbital compression discs and deleted components
        # correspond to either 3-spheres or edge-ideal triangulations of
        # unknots, so we can ignore all of the extra bookkeeping.
        decomposed, numOrbCuts, delComps, inconsistent = decomposeAlong(
                sphere, edgeIdealTri )
        newBlueprints = []
        for newEdgeIdealTri in decomposed:
            if isinstance( newEdgeIdealTri, EdgeIdealTriangulation ):
                try:
                    newEdgeIdealTri.simplify()
                except BoundsDisc:
                    # We can just throw away this unknot.
                    continue
                else:
                    # We are guaranteed to have len(newEdgeIdealTri) == 1.
                    newBlueprints.append( newEdgeIdealTri.blueprint() )
        sender.send(
                ( newBlueprints, attempts, searches, tri.size() ) )
        return
    return


def _enumerateSerial( oldEdgeIdealTri, tracker ):
    # Searching for quadrilateral vertex normal 2-spheres can be very slow.
    # However, if the oldEdgeIdealTri is a composite knot, then in practice
    # we find that we can "usually" find the desired 2-sphere very quickly.
    # Thus, when the enumeration takes a long time for the given
    # oldEdgeIdealTri, it is often helpful to randomise the edge-ideal
    # triangulation, and attempt the enumeration on the new triangulation.
    #
    # Unlike in _enumerateParallel(), here we implement the above idea in a
    # single-threaded fashion.
    tri = oldEdgeIdealTri.triangulation()
    enumeration = TreeEnumeration( tri, NS_QUAD )
    while True:
        if tracker.hasStalled():
            # We have spent a comparatively long time on the current
            # triangulation, so it might be worthwhile to try harder to
            # simplify this triangulation, and to restart the surface
            # enumeration on a smaller triangulation.
            tracker.report( None, "Try to simplify." )
            simpEdgeIdealTri = oldEdgeIdealTri.clone()
            simpEdgeIdealTri.randomise()    # Might raise BoundsDisc.
            if simpEdgeIdealTri.triangulation().size() < tri.size():
                oldEdgeIdealTri.setFromLoops(simpEdgeIdealTri)
                tri = oldEdgeIdealTri.triangulation()
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
            # No suitable 2-sphere means oldEdgeIdealTri represents a prime
            # knot.
            return None

        # We only want 2-spheres that intersect the oldEdgeIdealTri in either
        # exactly 0 points or exactly 2 points, since crushing such a
        # 2-sphere has one of the following effects:
        #   --> simplifies the triangulation containing the ideal loop;
        #   --> decomposes the oldEdgeIdealTri into two simpler knots; or
        #   --> (if oldEdgeIdealTri is unknotted) destroys all traces of the
        #       loop.
        wt = oldEdgeIdealTri.weight(sphere)
        if wt != 0 and wt != 2:
            continue
        # Because we are working with edge-ideal triangulations of knots in
        # the 3-sphere, orbital compression discs and deleted components
        # correspond to either 3-spheres or edge-ideal triangulations of
        # unknots, so we can ignore all of the extra bookkeeping.
        decomposed, numOrbCuts, delComps, inconsistent = decomposeAlong(
                sphere, EdgeIdealTriangulation( [oldEdgeIdealTri] ) )
        newEdgeIdealKnots = []
        for edgeIdealTri in decomposed:
            if isinstance( edgeIdealTri, EdgeIdealTriangulation ):
                try:
                    edgeIdealTri.simplify()
                except BoundsDisc:
                    # We can just throw away this unknot.
                    continue
                else:
                    # We are guaranteed to have len(edgeIdealTri) == 1.
                    newEdgeIdealKnots.append(edgeIdealTri)
        return newEdgeIdealKnots


def decomposeUsingAnnulus( knot, tracker=False, insertAsChild=False ):
    """
    Decomposes the given knot into prime pieces, represented as 3-spheres
    in which the prime knots are embedded as ideal loops.

    Unlike the decompose() routine, which works exclusively with ideal loops,
    the first step of this routine is to search for a quadrilateral vertex
    normal annulus in a triangulation of the knot exterior. If such an
    annulus exists, then crushing the annulus produces edge-ideal
    triangulations, and thereafter this routine also works entirely with
    ideal loops.

    The given knot should be provided as an instance of Regina's Link or
    PacketOfLink.

    If tracker is an instance of KnotDecompositionTracker, then this routine
    will use this given tracker to track the progress of the decomposition
    computation; if the tracker has the verbose option switched on, then this
    routine will also use the tracker to print regular progress reports.
    Otherwise, the routine will create its own KnotDecompositionTracker, and
    the tracker parameter should be either True or False depending on whether
    the newly-created tracker should have the verbose option switched on.

    If tracker is True or False, then the tracker created by this routine
    will have the timeout feature switched off. Thus, the only way to use the
    timeout feature with this routine is to explicitly supply a tracker.

    An explicitly supplied tracker may or may not already be started. If it
    is already started, then it will be assumed that the tracker is tracking
    a larger computation; hence, this routine will not call the tracker's
    finish() routine. Otherwise, if the tracker is not already started, then
    it will be assumed that the tracker is intended to track only the
    progress of this routine; hence, this routine will call the tracker's
    start() routine before performing the bulk of the computation, and it
    will also call the tracker's finish() routine after it has performed the
    bulk of the computation.

    If insertAsChild is True and the given knot is an instance of
    PacketOfLink, then this routine will insert the results of the
    computation as descendents of the given knot packet. This feature is also
    switched off by default.
    """
    #TODO Decide whether to bother with implementing decomposeUsingAnnulus().
    raise NotImplementedError()


class KnotDecompositionTracker:
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
        Creates a new KnotDecompositionTracker.

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

    def isStarted(self):
        """
        Has this tracker already been started?
        """
        return ( self._startTime is not None )

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

    def extendTimeout( self, seconds ):
        """
        Extends the allotted time by the given number of seconds.

        This routine does nothing if the timeout feature is switched off or
        the tracked computation has already finished.
        """
        if self._timeout is not None and self._finishTime is None:
            self._timeout += seconds
        return

    def checkTimeout(self):
        """
        Checks whether the tracked computation should be timed out, and if so
        raises TimeoutError.

        This routine does nothing if the timeout option is switched off or
        the allotted number of seconds has not yet elapsed.

        This routine must never be called before start() has been called.
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

    def newEdgeIdealTri( self, edgeIdealTri, extend=True ):
        """
        Informs this tracker that the tracked computation has started
        processing the given new edge-ideal triangulation.

        If this tracker is verbose, then this routine will automatically
        print a progress report.

        This routine raises TimeoutError if it detects that the tracked
        computation should be timed out.

        If extend is True (the default) and the timeout feature is switched
        on, then this routine extends the allotted time by a number of
        seconds equal to the number of tetrahedra in the given edge-ideal
        triangulation.

        This routine must never be called before start() has been called.
        """
        self._numTri += 1
        size = edgeIdealTri.triangulation().size()
        beforeReport = "Processing new {}-tetrahedron".format(size)
        beforeReport += " edge-ideal triangulation.\n"
        triEncoding, edgeIndices, orientation = edgeIdealTri.blueprint()
        beforeReport += "    Encoding:    {}\n".format(
                triEncoding.encode("unicode_escape") )
        beforeReport += "    Edges:       {}\n".format(
                ", ".join( [ str(i) for i in edgeIndices ] ) )
        beforeReport += "    Orientation: {}".format(orientation)

        # This counts as a new event.
        self._newEvent(beforeReport)
        if extend:
            self.extendTimeout(size)
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

    def unknownPrime( self, msg=None ):
        """
        Informs this tracker that the tracked computation has found a prime
        knot, but it is not yet known whether this prime is nontrivial.

        If this tracker is verbose, then this routine will automatically
        print a progress report. This progress report may be preceded by an
        additional message supplied via the optional msg argument.

        This routine raises TimeoutError if it detects that the tracked
        computation should be timed out.

        This routine must never be called before start() has been called.
        """
        afterReport = "Found a prime knot! Is it nontrivial?"
        if msg is not None:
            afterReport = "{}\n{}".format( msg, afterReport )
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
