"""
Decompose knots into prime knots.
"""
from sys import stdout
from timeit import default_timer
from regina import *
from idealedge import decomposeAlong, isSphere
from loop import IdealLoop, BoundsDisc
from knotted import isKnotted, knownHyperbolic
from embed import loopPacket, reversePinch, embedByFilling, embedFromDiagram
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
    if isinstance( tracker, DecompositionTracker ):
        verbose = tracker.isVerbose()
    else:
        verbose = bool(tracker)
        tracker = DecompositionTracker(verbose)
    if tracker.isStarted():
        needToFinish = False
    else:
        needToFinish = True
        tracker.start()

    # Build the IdealLoop on which we perform the decomposition computation.
    # Make sure to create clones so as not to directly modify the input.
    if isinstance( knot, IdealLoop ):
        loop = knot.clone()
    elif isinstance( knot, Edge3 ):
        loop = IdealLoop( [knot] ).clone()
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
                loop = embedFromDiagram(knot)
            except BoundsDisc:
                # The given knot is unknotted.
                loop = None
        else:
            try:
                loop = _embedParallel( knot, tracker )
            except BoundsDisc:
                # The given knot is unknotted.
                loop = None

    # Do the decompositon.
    primes = []
    if loop is None:
        # The given knot is unknotted.
        if verbose:
            afterReport = "The knot bounds a disc!"
            tracker.report( None, afterReport )
        toProcess = []
    else:
        toProcess = [loop]
    while toProcess:
        # INVARIANT:
        #   At this point, the following are guaranteed to hold:
        #   --> Each element of toProcess is an ideal loop forming a
        #       (possibly trivial, possibly composite) knot.
        #   --> Each element of primes is an ideal loop forming a nontrivial
        #       prime knot.
        #   --> The input knot is given by composing all of the knots
        #       represented in toProcess and primes.
        oldLoop = toProcess.pop()
        tracker.newLoop(oldLoop)
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
            try:
                newLoops = _enumerateSerial( oldLoop, tracker )
            except BoundsDisc:
                # The oldLoop is unknotted.
                tracker.knownPrime(False)
                continue
            msg = None
        else:
            try:
                newLoops, msg = _enumerateParallel( oldLoop, tracker )
            except BoundsDisc:
                # The oldLoop is unknotted.
                tracker.knownPrime(False)
                continue
        if newLoops is None:
            # The oldLoop is prime! However, we only care about the case
            # where this prime is nontrivial.
            tracker.unknownPrime(msg)
            isNontrivial = isKnotted( oldLoop, tracker )
            if isNontrivial:
                primes.append(oldLoop)
            tracker.knownPrime(isNontrivial)
        else:
            toProcess.extend(newLoops)
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
        for i, primeLoop in enumerate(primes):
            packet = loopPacket(primeLoop)
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

        # Have we finished embedding the knot as an ideal loop?
        if not fillingProcess.is_alive():
            diagramProcess.terminate()
            fillingProcess.join()
            diagramProcess.join()
            if fillingReceiver.poll():
                loop = IdealLoop()
                loop.setFromLightweight( *fillingReceiver.recv() )
            else:
                # If fillingProcess terminated without sending information,
                # then the given knot must be unknotted.
                raise BoundsDisc()
            if tracker is not None:
                afterReport = "Built triangulation using 1/0 Dehn surgery."
                tracker.report( None, afterReport )
            return loop
        if not diagramProcess.is_alive():
            fillingProcess.terminate()
            diagramProcess.join()
            fillingProcess.join()
            if diagramReceiver.poll():
                loop = IdealLoop()
                loop.setFromLightweight( *diagramReceiver.recv() )
            else:
                # If diagramProcess terminated without sending information,
                # then the given knot must be unknotted.
                raise BoundsDisc()
            if tracker is not None:
                afterReport = "Built triangulation from planar diagram."
                tracker.report( None, afterReport )
            return loop
    return


def _runFilling( knotSig, sender ):
    RandomEngine.reseedWithHardware()
    try:
        loop = embedByFilling( Link.fromKnotSig(knotSig) )
    except BoundsDisc:
        # Send nothing if the given knot is unknotted.
        return
    sender.send( loop.lightweightDescription() )
    return


def _runDiagram( knotSig, sender ):
    RandomEngine.reseedWithHardware()
    try:
        loop = embedFromDiagram( Link.fromKnotSig(knotSig) )
    except BoundsDisc:
        # Send nothing if the given knot is unknotted.
        return
    sender.send( loop.lightweightDescription() )
    return


def _enumerateParallel( oldLoop, tracker ):
    # Searching for quadrilateral vertex normal 2-spheres can be very slow.
    # However, if the oldLoop is a composite knot, then in practice we find
    # that we can "usually" find the desired 2-sphere very quickly. Thus,
    # when the enumeration takes a long time for the given oldLoop, it is
    # often helpful to randomise the loop and attempt the enumeration on the
    # new loop.
    description = oldLoop.lightweightDescription()
    tri = oldLoop.triangulation()

    # Set up a child process to repeatedly randomise the given ideal loop,
    # and send the randomised loops to another child process that runs
    # alternate enumerations.
    randomiseReceiver, randomiseSender = Pipe(False)
    randomiseProcess = Process( target=_perpetualRandomise,
            args=( description, tri.size(), randomiseSender ) )
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
        # Has the randomiseProcess determined that the oldLoop is unknotted?
        if not randomiseProcess.is_alive():
            # Make sure to clean up child processes before raising BoundsDisc
            # to indicate that the oldLoop is unknotted.
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
            newLoopDescs, attempts, searches, size = alternateReceiver.recv()
            msg = "Alternate enumeration succeeded on "
            msg += "{}-tetrahedron triangulation.\n".format(size)
            msg += "(Randomisation attempts: {}. Searches: {}.)".format(
                    attempts, searches )
            if newLoopDescs is None:
                # Found a prime!
                return ( None, msg )
            else:
                # Build new loops and return them.
                newLoops = []
                for description in newLoopDescs:
                    newLoop = IdealLoop()
                    newLoop.setFromLightweight( *description )
                    newLoops.append(newLoop)
                return ( newLoops, msg )

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
            # No suitable 2-sphere means oldLoop is prime.
            # Clean up child processes before returning.
            alternateProcess.terminate()
            randomiseProcess.terminate()
            alternateProcess.join()
            randomiseProcess.join()
            return ( None, msg )

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

        # Clean up child processes before returning.
        alternateProcess.terminate()
        randomiseProcess.terminate()
        alternateProcess.join()
        randomiseProcess.join()
        return ( newLoops, msg )
    return


def _perpetualRandomise( description, size, sender ):
    RandomEngine.reseedWithHardware()
    loop = IdealLoop()
    loop.setFromLightweight( *description )
    attempts = 0
    while True:
        attempts += 1
        try:
            loop.randomise()    # Might raise BoundsDisc.
        except BoundsDisc:
            # Use early termination to indicate that the loop is unknotted.
            return
        if loop.triangulation().size() <= size:
            # Send randomised loop.
            sender.send( ( loop.lightweightDescription(), attempts ) )
    return


def _indefiniteEnumerate( receiver, sender ):
    loop = IdealLoop()
    searches = 0
    while not receiver.poll():
        sleep(0.01)
    description, attempts = receiver.recv()
    loop.setFromLightweight( *description )
    tri = loop.triangulation()
    enumeration = TreeEnumeration( tri, NS_QUAD )
    while True:
        if searches > 20 and receiver.poll():
            # Restart the enumeration with a new ideal loop.
            searches = 0
            description, attempts = receiver.recv()
            loop.setFromLightweight( *description )
            tri = loop.triangulation()
            enumeration = TreeEnumeration( tri, NS_QUAD )

        # Get the next 2-sphere.
        searches += 1
        if enumeration.next():
            sphere = enumeration.buildSurface()
            if not isSphere(sphere):
                continue
        else:
            # No suitable 2-sphere means the loop is prime.
            sender.send( ( None, attempts, searches, tri.size() ) )
            return

        # We only want 2-spheres that intersect the loop in either exactly 0
        # points or exactly 2 points, since crushing such a 2-sphere has one
        # of the following effects:
        #   --> simplifies the triangulation containing the ideal loop;
        #   --> decomposes the loop into two simpler knots; or
        #   --> (if the loop is unknotted) destroys all traces of the loop.
        wt = loop.weight(sphere)
        if wt != 0 and wt != 2:
            continue
        decomposed = decomposeAlong( sphere, [loop] )
        newLoopDescriptions = []
        for decomposedLoops in decomposed:
            if decomposedLoops:
                # We are guaranteed to have len(decomposedLoops) == 1.
                newLoopDescriptions.append(
                        decomposedLoops[0].lightweightDescription() )
        sender.send(
                ( newLoopDescriptions, attempts, searches, tri.size() ) )
        return
    return


def _enumerateSerial( oldLoop, tracker ):
    # Searching for quadrilateral vertex normal 2-spheres can be very slow.
    # However, if the oldLoop is a composite knot, then in practice we find
    # that we can "usually" find the desired 2-sphere very quickly. Thus,
    # when the enumeration takes a long time for the given oldLoop, it is
    # often helpful to randomise the loop and attempt the enumeration on the
    # new loop.
    #
    # Unlike in _enumerateParallel(), here we implement the above idea in a
    # single-threaded fashion.
    tri = oldLoop.triangulation()
    enumeration = TreeEnumeration( tri, NS_QUAD )
    while True:
        if tracker.hasStalled():
            # We have spent a comparatively long time on the current
            # triangulation, so it might be worthwhile to try harder to
            # simplify this triangulation, and to restart the surface
            # enumeration on a smaller triangulation.
            tracker.report( None, "Try to simplify." )
            simpLoop = oldLoop.clone()
            simpLoop.randomise()    # Might raise BoundsDisc.
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
            return None

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
        return newLoops


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
    #TODO
    raise NotImplementedError()


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

    def newLoop( self, loop, extend=True ):
        """
        Informs this tracker that the tracked computation has started
        processing the given new ideal loop.

        If this tracker is verbose, then this routine will automatically
        print a progress report.

        This routine raises TimeoutError if it detects that the tracked
        computation should be timed out.

        If extend is True (the default) and the timeout feature is switched
        on, then this routine extends the allotted time by a number of
        seconds equal to the size of the triangulation containing the given
        ideal loop.

        This routine must never be called before start() has been called.
        """
        self._numTri += 1
        size = loop.triangulation().size()
        sig, edgeLocations = loop.lightweightDescription()
        beforeReport = "Edge-ideal: {} ".format(sig)
        if size == 1:
            beforeReport += "(1 tetrahedron). "
        else:
            beforeReport += "({} tetrahedra). ".format(size)

        # Work out edge indices after reconstructing from iso sig.
        temp = Triangulation3.fromIsoSig(sig)
        edgeIndices = []
        for tetIndex, edgeNum in edgeLocations:
            edgeIndices.append( str(
                    temp.tetrahedron(tetIndex).edge(edgeNum).index() ) )
        if len(loop) == 1:
            beforeReport += "Edge {}.".format( edgeIndices[0] )
        else:
            beforeReport += "Edges {}.".format( ", ".join(edgeIndices) )

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
