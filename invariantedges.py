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

    def _setFromEdgeLocationsImpl( self, edgeLocations ):
        """
        Sets this embedded loop using the given collection of edge locations.

        In detail, each edge location must be a pair (t, e), where t is a
        tetrahedron and e is an edge number from 0 to 5 (inclusive). Each
        tetrahedron must belong to the same triangulation.

        This routine is *not* implemented by default, so subclasses that
        require this routine must provide an implementation. Subclasses are
        free to specify the data structure for the collection of edge
        locations, as long as each edge location is represented as a pair of
        the form described above.
        """
        raise NotImplementedError()

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
        """
        Shortens this collection of invariant edges by replacing it with a
        smaller collection of edges that is topologically equivalent to the
        current collection.

        The default implementation of this routine repeatedly searches for a
        triangle F that intersects exactly two of the edges in this
        collection, and redirects these two edges along the third edge of F.
        Only the logic of the search-and-redirect loop is explicitly
        implemented, and the rest is deferred to the following helper
        routines, which are *not* implemented by default:
        --> _findRedirect(), which finds an available redirect operation (if
            any); and
        --> _performRedirect(), which performs a redirect operation.

        Subclasses that require this routine must therefore either:
        --> override this routine; or
        --> supply implementatiions for the helper routines listed above.
        In the latter case, see the documentation for these helper routines
        for details on the behaviour that must be implemented.

        Any implementation of this routine must return True if shortening was
        successful, and False otherwise. Moreover, in the False case, the
        collection of invariant edges must be left entirely untouched.
        """
        changed = False     # Eventual return value.

        # Continue attempting to redirect until no further such operations
        # are possible.
        while True:
            redirect = self._findRedirect()
            if redirect is None:
                break
            self._performRedirect(*redirect)
            changed = True
        return changed

    def _findRedirect(self):
        """
        Finds an edge of a triangle that can be used as the target of a
        redirect operation, or returns None if no such redirect operations
        can be performed.

        In detail, consider an edge e of a triangle F. If the other two edges
        of F are distinct and are both part of this collection of invariant
        edges, but e is not part of this collection, then the redirect
        operation consists of removing the other two edges from this
        collection and replacing them with e. This operation can be legally
        performed if the new collection of invariant edges is topologically
        equivalent to the current collection. In the case where this
        operation is indeed legal, this routine returns the pair (F, n),
        where n is the edge number of F corresponding to the edge e.

        The InvariantEdges base class does not implement this routine, so
        subclasses that require this routine must provide an implementation.
        """
        raise NotImplementedError()

    def _performRedirect( self, triangle, edgeNum ):
        """
        Performs a redirect across the given triangle that targets the edge
        triangle.edge(edgeNum).

        Pre-condition:
        --> The requested redirect operation must be legal, as described in
            the documentation for the _findRedirect() routine.
        """
        raise NotImplementedError()

    def _minimiseBoundaryImpl(self):
        """
        Attempts to minimise the number of boundary triangles in the ambient
        triangulation without topologically altering this collection of
        invariant edges.

        The default implementation of this routine attempts to reduce the
        number of boundary triangles by performing the following operations
        until no further such operations are possible:
        --> Shortening this collection of invariant edges using the
            _shortenImpl() routine.
        --> Performing close book moves on the ambient triangulation.
        --> Layering across a boundary edge of the ambient triangulation, and
            then immediately performing a close book move across the new
            boundary edge obtained after layering.
        Much of this implementation is deferred to the following helper
        routines, which are *not* fully implemented by default:
        --> _shortenImpl()
        --> _findBoundaryMove()
        --> _setFromEdgeLocationsImpl()

        Subclasses that require this routine must therefore either:
        --> override this routine; or
        --> supply implementatiions for the helper routines listed above.
        In the latter case, see the documentation for these helper routines
        for details on the behaviour that must be implemented.

        Adapted from Regina's Triangulation3.minimiseBoundary().

        Precondition:
        --> The ambient triangulation is valid.

        Returns:
            True if and only if a change was made to either this collection
            of invariant edges or its ambient triangulation (or both).
        """
        changed = False     # Eventual return value.
        while True:
            if self._shortenImpl():
                changed = True

            # Is there a move we can perform to reduce the number of boundary
            # triangles?
            boundaryMove = self._findBoundaryMove()
            if boundaryMove is None:
                # At this point, there is nothing "easy" we can do to further
                # reduce the number of boundary triangles.
                return changed

            # Perform the move, and then update this collection of invariant
            # edges. The _findBoundaryMove() routine should have already
            # checked that the close book move is legal.
            changed = True
            edge, doLayer, newEdgeLocations = boundaryMove
            if doLayer:
                edge = self._tri.layerOn(edge).edge(5)
            self._tri.closeBook( edge, False, True )
            self._setFromEdgeLocationsImpl(edgeLocations)

        # Should never reach this point.
        raise RuntimeError(
                "Unexpectedly broke out of boundary move loop." )

    def _findBoundaryMove(self):
        """
        Returns details of a boundary move that reduces the number of
        boundary triangles in the ambient triangulation, or None if no such
        moves are available.

        In detail, this routine checks whether (at least) one of the
        following operations is possible, subject to the constraint that this
        collection of invariant edges is not topologically altered:
        --> Performing close book moves on the ambient triangulation.
        --> Layering across a boundary edge of the ambient triangulation, and
            then immediately performing a close book move across the new
            boundary edge obtained after layering.

        When such a boundary move *does* exist, the return value will be a
        tuple that describes the move using the following data:
        (0) A boundary edge e on which to perform the move.
        (1) A boolean indicating whether we need to layer across e before
            performing the close book move.
        (2) A collection of edge locations that describes the location of the
            invariant edges after performing the move. More precisely, after
            the move, the collection of invariant edges can be found by
            passing the edge locations to _setFromEdgeLocationsImpl().

        The InvariantEdges base class does not implement this routine, so
        subclasses that require this routine must provide an implementation.
        """
        raise NotImplementedError()

    #TODO

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
            self._setFromEdgeLocationsImpl(edgeLocations)
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
