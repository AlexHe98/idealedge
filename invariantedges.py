"""
A class to track edges that should be kept invariant as we simplify the
ambient triangulation.
"""
from regina import *
from moves import twoThree, threeTwo, twoZero, twoOne, fourFour


def snapEdge( edge, check=True, perform=True ):
    """
    If the endpoints of the given edge are distinct and not both boundary,
    then uses a snapped ball to pinch these two endpoints together.

    This operation is equivalent to performing the following two operations:
    (1) Pinching the edge, which introduces a two-tetrahedron gadget with a
        single degree-one edge e at its heart.
    (2) Performing a 2-1 edge move on e.

    If check is True (the default), then this routine will check whether
    snapping the given edge is legal; otherwise, this routine will proceed
    under the assumption that the move is already known to be legal. If
    perform is True (the default), then this routine will actually perform
    the snap edge move if it has determined or assumed that the move is
    legal; otherwise, the triangulation containing the given edge will be
    left unchanged.

    If the triangulation containing the given edge is currently oriented,
    then this operation will preserve the orientation.

    Parameters:
    --> edge    The edge whose endpoints should be snapped together.

    Returns:
        True if and only if snapping the given edge is legal.
    """
    if check:
        # Endpoints need to be distinct and not both boundary.
        u = edge.vertex(0)
        v = edge.vertex(1)
        if u == v:
            return False
        if u.isBoundary() and v.isBoundary():
            return False
    if not perform:
        return True

    # Start by pinching the given edge.
    tri = edge.triangulation()
    tri.pinchEdge(edge)

    # To find the degree-one edge at the heart of the pinch edge gadget, look
    # at the last two tetrahedra in tri.
    found = False
    for tetIndex in [ tri.size() - 1, tri.size() - 2 ]:
        for edgeNum in range(6):
            e = tri.tetrahedron(tetIndex).edge(edgeNum)
            if e.degree() == 1:
                found = True
                break
        if found:
            break

    # Finish up by performing a 2-1 move on e.
    if not tri.twoOneMove( e, 0 ):
        if not tri.twoOneMove( e, 1 ):
            # This should never happen.
            raise RuntimeError( "Snap edge failed unexpectedly." )
    return True


class InvariantEdges:
    """
    A collection of edges that should be kept invariant as we simplify the
    ambient triangulation.
    """
    def __init__(self):
        self._tri = None
        self._edgeIndices = None
        return

    def _setTriangulation( self, tri ):
        """
        Sets the ambient triangulation.
        """
        self._tri = tri
        return

    def _setEdgeIndices( self, edgeIndices ):
        """
        Sets the indices of the invariant edges.
        """
        self._edgeIndices = edgeIndices
        return

    def triangulation(self):
        """
        Returns the ambient triangulation.
        """
        return self._tri

    def __len__(self):
        return len( self._edgeIndices )

    def __contains__( self, edgeIndex ):
        return edgeIndex in self._edgeIndices

    def __iter__(self):
        return iter( self._edgeIndices )

    def intersects( self, surf ):
        """
        Returns True if and only if this collection of invariant edges has
        nonempty intersection with the given normal surface surf.

        Precondition:
        --> The given normal surface is embedded in self.triangulation().
        """
        for i in self:
            if surf.edgeWeight(i).safeLongValue() > 0:
                return True
        return False

    def weight( self, surf ):
        """
        Returns the number of times this collection of invariant edges
        intersects the given normal surface surf.

        Precondition:
        --> The given normal surface is embedded in self.triangulation().
        """
        wt = 0
        for i in self:
            wt += surf.edgeWeight(i).safeLongValue()
        return wt

    def _shortenImpl(self):
        #TODO
        """
        Shortens this collection of invariant edges by replacing it with a
        smaller collection of edges that is topologically equivalent to the
        current collection.

        Topological equivalence depends on context, so this routine is not
        implemented by default. Thus, subclasses that require this routine
        must supply their own implementations.

        Implementations must return True if shortening was successful, and
        False otherwise. Moreover, in the False case, the collection of
        invariant edges must be left entirely untouched.
        """
        raise NotImplementedError()

    def _minimiseBoundaryImpl(self):
        #TODO
        """
        Ensures that the triangulation containing this collection of
        invariant edges has the smallest possible number of boundary
        triangles, potentially adding tetrahedra to do this.

        The default implementation of this routine requires the following
        helper routines, which are *not* implemented by default:
        --> _shortenImpl()
        --> _findBoundaryMove()
        Thus, subclasses that require this routine must either:
        --> override this routine; or
        --> supply implementations for the aforementioned helper routines.
        In the latter case, see the documentation for each respective helper
        routine for details on the behaviour that must be implemented.

        The changes that this routine performs can always be expressed using
        only the following operations:
        --> Shortening by running the _shortenImpl() routine.
        --> Close book moves and/or layerings on self.triangulation().
        In particular, this routine never creates new vertices, and it never
        creates a non-vertex-linking normal disc or 2-sphere if there was not
        one before.

        After this routine is finished, it is guaranteed that, subject to the
        constraint that this collection of invariant edges must be preserved,
        it will be impossible to simplify the boundary of self.triangulation()
        any further using the operations listed above.

        Adapted from Regina's Triangulation3.minimiseBoundary().

        Precondition:
        --> The ambient triangulation (i.e., self.triangulation()) is valid.

        Returns:
            True if and only if this loop or its ambient triangulation were
            changed. In other words, a return value of False indicates that:
            (1) this collection of invariant edges could not be shortened;
                and
            (2) every boundary component of the ambient triangulation was
                already minimal to begin with.
        """
        #TODO
        # Simplify the boundary of self._tri to use as few triangles as
        # possible.
        changed = False
        while True:
            # Shorten to minimise the number of special cases.
            if self._shortenImpl():
                changed = True

            # Is there a move we can perform to reduce the number of boundary
            # triangles?
            moveDetails = self._findBoundaryMove()
            if moveDetails is None:
                # At this point, all boundary components must be minimal.
                return changed
            changed = True
            edge, doLayer, newEdgeIndices = moveDetails

            # Make sure we will be able to find the edges that form the loop
            # after performing the move.
            edgeLocations = []
            for ei in newEdgeIndices:
                emb = self._tri.edge(ei).embedding(0)
                edgeLocations.append( ( emb.tetrahedron(), emb.edge() ) )

            # Perform the move, and then update this loop. We can safely
            # assume that the close book move is legal, so no need to check
            # before performing.
            if doLayer:
                edge = self._tri.layerOn(edge).edge(5)
            self._tri.closeBook( edge, False, True )
            self.setFromEdgeLocations(edgeLocations)
        return

    #TODO

    def _findBoundaryMove(self):
        """
        Returns details of a boundary move that simplifies the boundary of
        self.triangulation(), or None if the boundary is already minimal.

        In detail, in the case where the boundary is not yet minimal, this
        routine guarantees to find a move that reduces the number of boundary
        triangles by two (without changing the topology of this loop). The
        return value will be a tuple that describes this move using the
        following data:
        (0) A boundary edge e on which to perform the move.
        (1) A boolean indicating whether we need to layer across e. If this
            is True, then the move we perform will be to first layer across
            e, and then perform a close book move on the newly layered edge.
            Otherwise, the move will simply be a close book move on e.
        (2) A list of edge indices describing a sequence of edges (in the
            current triangulation) such that, after performing the move, this
            edge sequence becomes a loop that is topologically equivalent to
            this loop.

        The EmbeddedLoop base class does not implement this routine, so
        subclasses that require this routine must provide an implementation.
        """
        raise NotImplementedError()

    def _minimiseVerticesImpl(self):
        """
        Ensures that the triangulation containing this embedded loop has the
        smallest possible number of vertices for the 3-manifold that it
        represents, potentially adding tetrahedra to do this.

        The default implementation of this routine requires the following
        helper routines, which are *not* fully implemented by default:
        --> _shortenImpl()
        --> _minimiseBoundaryImpl()
        --> _findSnapEdge()
        Thus, subclasses that require this routine must either:
        --> override this routine; or
        --> supply implementations for the aforementioned helper routines.
        In the latter case, see the documentation for each respective helper
        routine for details on the behaviour that must be implemented.

        A side-effect of calling this routine is that it will shorten this
        embedded loop if possible.

        This routine might raise BoundsDisc.

        The following are guaranteed to hold once this routine is finished:
        --> If the ambient triangulation is closed, then it will have
            precisely one vertex.
        --> If the ambient triangulation has real boundary, then:
            --- either there will be no internal vertices, or there will be
                exactly one internal vertex if this embedded loop is required
                to remain entirely in the interior of the triangulation;
            --- every 2-sphere boundary component will have exactly two
                triangles and three vertices;
            --- every projective plane boundary component will have exactly
                two triangles and two vertices;
            --- every other boundary component will have exactly one vertex.

        The changes that this routine performs can always be expressed using
        only the following operations:
        --> Shortening this loop by redirecting it across triangular faces.
        --> Close book moves, layerings and/or snap edge moves on
            self.triangulation().
        In particular, this routine never creates new vertices.

        Adapted from Regina's Triangulation3.minimiseVertices().

        Precondition:
        --> The ambient triangulation (i.e., self.triangulation()) is valid.

        Returns:
            True if and only if this loop or its ambient triangulation were
            changed. In other words, a return value of False indicates that:
            (1) this loop could not be shortened; and
            (2) the number of vertices in the ambient triangulation was
                already minimal to begin with.
        """
        # Start by minimising the boundary.
        changed = self._minimiseBoundaryImpl()  # Might raise BoundsDisc.

        # All that remains now is to remove internal vertices.
        # We do not currently have an implementation of collapseEdge() that
        # keeps track of how edges get relabelled after the move, so we rely
        # entirely on the snapEdge() routine.
        while True:
            # Shorten this loop to minimise the number of special cases.
            if self._shortenImpl():     # Might raise BoundsDisc.
                changed = True

            # Is there a snap edge move we can perform to reduce the number
            # of vertices?
            moveDetails = self._findSnapEdge()
            if moveDetails is None:
                # At this point, there are no more unnecessary internal
                # vertices.
                return changed
            changed = True
            edge, newEdgeIndices = moveDetails

            # Make sure we will be able to find the edges that form the loop
            # after performing the move.
            edgeLocations = []
            for ei in newEdgeIndices:
                emb = self._tri.edge(ei).embedding(0)
                edgeLocations.append( ( emb.tetrahedron(), emb.edge() ) )

            # Perform the snap, and then update this ideal loop. We can
            # assume that the snap is legal, so can perform without checking.
            snapEdge( edge, False, True )
            self.setFromEdgeLocations(edgeLocations)
        return

    def _findSnapEdge(self):
        """
        Returns details of a snap edge move that can be used to reduce the
        number of vertices() in self.triangulation(), or None if the number
        of vertices is already minimal.

        In detail, in the case where the number of vertices is not yet
        minimal, this routine returns a tuple consisting of the following:
        (0) An edge on which a snap edge move can be performed.
        (1) A list of edge indices describing a sequence of edges (in the
            current triangulation) such that, after performing the move, this
            edge sequence becomes a loop that is topologically equivalent to
            this loop.

        The EmbeddedLoop base class does not implement this routine, so
        subclasses that require this routine must provide an implementation.

        Pre-condition:
        --> This loop cannot be shortened.
        --> If the ambient triangulation has real boundary, then this
            boundary has already been minimised.
        """
        raise NotImplementedError()

    def _simplifyMonotonicImpl( self, include32 ):
        """
        Uses 2-0 edge, 2-1 edge, and (optionally) 3-2 moves to monotonically
        reduce the number of tetrahedra in the ambient triangulation, while
        leaving this embedded loop untouched.

        This routine might raise BoundsDisc.

        Adapted from Regina's Triangulation3.simplifyToLocalMinimum() and
        SnapPea's check_for_cancellation().

        Returns:
            True if and only if the ambient triangulation was successfully
            simplified. Otherwise, the ambient triangulation will not be
            modified at all.
        """
        changed = False     # Has anything changed ever?    (Return value.)
        changedNow = True   # Did we just change something? (Loop control.)
        while True:
            changedNow = False
            for edge in self._tri.edges():
                # Make sure to leave the embedded loop untouched.
                if edge.index() in self:
                    continue

                # If requested, try a 3-2 move.
                if include32:
                    renum = threeTwo(edge)
                    if renum is not None:
                        changedNow = True
                        changed = True
                        break

                # Try a 2-0 edge move.
                # This move can destroy the loop if it bounds a disc.
                renum = twoZero(edge)
                if renum is not None:
                    changedNow = True
                    changed = True
                    break

                # Try a 2-1 edge move.
                # This move can destroy the loop if it bounds a disc.
                renum = twoOne( edge, 0 )
                if renum is not None:
                    changedNow = True
                    changed = True
                    break
                renum = twoOne( edge, 1 )
                if renum is not None:
                    changedNow = True
                    changed = True
                    break

            # Did we improve the triangulation? If so, then we need to update
            # the details of the loop, and then check whether we can make
            # further improvements.
            if changedNow:
                try:
                    # If we destroyed the loop, then this will raise NotLoop.
                    self._setFromRenum(renum)
                except NotLoop:
                    # As noted above, the loop can only get destroyed if it
                    # bounds a disc.
                    raise BoundsDisc()
            else:
                break

        # Nothing further we can do.
        return changed

    def _simplifyImpl(self):
        """
        Attempts to simplify this embedded loop.

        This routine uses _minimiseVerticesImpl() and
        _simplifyMonotonicImpl(), in combination with random 4-4 moves that
        leave this loop untouched.

        Note that the helper routine _minimiseVerticesImpl() is *not* fully
        implemented by default. Thus, subclasses that require this
        _simplifyImpl() routine must either:
        --> override this routine; or
        --> supply an implementation for _minimiseVerticesImpl().
        In the latter case, see the documentation for _minimiseVerticesImpl()
        for details on the behaviour that must be implemented.

        This routine might raise BoundsDisc.

        Adapted from Regina's Triangulation3.intelligentSimplify().

        Warning:
            --> Running this routine multiple times on the same loop may
                return different results, since the implementation makes
                random decisions.

        Returns:
            True if and only if this loop was successfully simplified.
            Otherwise, this loop will not be modified at all.
        """
        RandomEngine.reseedWithHardware()

        # Work with a clone so that we can roll back changes if we are not
        # able to reduce the number of tetrahedra.
        tempLoop = self.clone()
        tempTri = tempLoop.triangulation()

        # Start by minimising vertices. This will probably increase the
        # number of tetrahedra if the number of vertices is not already
        # minimal, but hopefully the monotonic simplification saves us.
        #
        # Might raise BoundsDisc.
        tempLoop._minimiseVerticesImpl()
        tempLoop._simplifyMonotonicImpl(True)   # Include 3-2 moves.

        # Use random 4-4 moves until it feels like even this is not helping
        # us make any further progress.
        #
        # In detail, we keep track of a cap on the number of consecutive 4-4
        # moves that we are allowed to perform without successfully
        # simplifying the triangulation. We give up whenever we reach or
        # exceed this cap. The cap is scaled up based on our "perseverance".
        fourFourAttempts = 0
        fourFourCap = 0
        perseverance = 5        # Hard-coded value copied from Regina.
        while True:
            # Find all available 4-4 moves.
            fourFourAvailable = []
            for edge in tempTri.edges():
                if edge.index() in tempLoop:
                    # We do not want to touch the embedded loop.
                    continue
                for axis in range(2):
                    if tempTri.fourFourMove( edge, axis, True, False ):
                        fourFourAvailable.append( ( edge, axis ) )

            # Is it worthwhile to continue attempting 4-4 moves?
            availableCount = len(fourFourAvailable)
            if fourFourCap < perseverance * availableCount:
                fourFourCap = perseverance * availableCount
            if fourFourAttempts >= fourFourCap:
                break

            # Perform a random 4-4 move, and see if this is enough to help us
            # simplify the triangulation.
            #
            # _simplifyMonotonicImpl() might raise BoundsDisc.
            fourFourChoice = fourFourAvailable[
                    RandomEngine.rand(availableCount) ]
            renum = fourFour( *fourFourChoice )
            tempLoop._setFromRenum(renum)
            if tempLoop._simplifyMonotonicImpl(True):   # Include 3-2 moves.
                # We successfully simplified!
                # Start all over again.
                fourFourAttempts = 0
                fourFourCap = 0
            else:
                fourFourAttempts += 1

        # If simplification was successful (either by reducing the number of
        # tetrahedra, or failing that by reducing the number of vertices
        # without increasing the number of tetrahedra), then sync this
        # embedded loop with the now-simpler tempLoop.
        tri = self.triangulation()
        simplified = ( tempTri.size() < tri.size() or
                ( tempTri.size() == tri.size() and
                    tempTri.countVertices() < tri.countVertices() ) )
        if simplified:
            self.setFromLoop(tempLoop)
        return simplified