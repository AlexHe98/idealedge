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
    A collection of edges that should be kept invariant (i.e., we must
    preserve the topological embedding of this collection of edges) as we
    simplify the ambient triangulation.

    This is a base class that provides some basic implementations that are
    intended to be inherited and/or extended by more specialised types of
    collections of edges.
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

    def setAsCloneOf( self, other, copyTri=True ):
        """
        Sets this collection of invariant edges to be a clone of the other
        collection of invariant edges.

        If copyTri is True (the default), then this collection of invariant
        edges will be embedded in a copy of the triangulation containing the
        other collection. Otherwise, if copyTri is False, then this
        collection of invariant edges will be embedded in the same
        triangulation as the other collection.

        The InvariantEdges base class does not implement this routine, so
        subclasses that require this routine must provide an implementation.
        """
        raise NotImplementedError()

    def _setFromEdgeLocationsImpl( self, edgeLocations ):
        """
        Sets this embedded loop using the given collection of edge locations.

        In detail, each edge location must be a pair (t, e), where t is a
        tetrahedron and e is an edge number from 0 to 5 (inclusive). Each
        tetrahedron must belong to the same triangulation.

        The InvariantEdges base class does not implement this routine, so
        subclasses that require this routine must provide an implementation.

        Such an implementation is allowed (but is not required) to include
        the following features:
        --> It may specify any particular data structure for the collection
            of edge locations, as long as each edge location is represented
            as a pair of the form described above.
        --> It may raise an error if it detects that the given edge locations
            have some property that is undesirable (what "undesirable" means
            may depend on the particular requirements of the subclass).
        """
        raise NotImplementedError()

    def clone(self):
        """
        Returns a clone of this collection of invariant edges.

        The cloned collection of invariant edges will always be embedded in a
        copy of the triangulation containing this collection of invariant
        edges.

        The InvariantEdges base class does not implement this routine, so
        subclasses that require this routine must provide an implementation.
        """
        raise NotImplementedError()

    def triangulation(self):
        """
        Returns the ambient triangulation.
        """
        return self._tri

    def __len__(self):
        """
        How many edges are there in this collection of invariant edges?

        The InvariantEdges base class does not implement this routine, so
        subclasses that require this routine must provide an implementation.
        """
        raise NotImplementedError()

    def __contains__( self, edgeIndex ):
        """
        Does this collection of invariant edges contain the edge with the
        given index?

        The InvariantEdges base class does not implement this routine, so
        subclasses that require this routine must provide an implementation.
        """
        raise NotImplementedError()

    def __iter__(self):
        """
        Returns an iterator that runs through all the edges in this
        collection of invariant edges.

        The InvariantEdges base class does not implement this routine, so
        subclasses that require this routine must provide an implementation.
        """
        raise NotImplementedError()

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
        implemented, and the rest is deferred to the following subroutines,
        which are *not* implemented by default:
        --> _findRedirect(), which finds an available redirect operation (if
            any); and
        --> _performRedirect(), which performs a redirect operation.

        Subclasses that require this routine must therefore either:
        --> override this routine; or
        --> supply implementations for the subroutines listed above.
        In the latter case, see the documentation for these subroutines for
        details on the behaviour that must be implemented.

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

        The InvariantEdges base class does not implement this routine, so
        subclasses that require this routine must provide an implementation.

        Pre-condition:
        --> The requested redirect operation must be legal, as described in
            the documentation for the _findRedirect() routine.
        """
        raise NotImplementedError()

    def _minimiseBoundaryImpl(self):
        """
        Attempts to minimise the number of boundary triangles in the ambient
        triangulation without topologically altering this collection of
        invariant edges, potentially adding tetrahedra to do this.

        The default implementation of this routine attempts to reduce the
        number of boundary triangles by performing the following operations
        until no further such operations are possible:
        --> Shortening this collection of invariant edges using the
            _shortenImpl() routine.
        --> Performing close book moves on the ambient triangulation.
        --> Layering across a boundary edge of the ambient triangulation, and
            then immediately performing a close book move across the new
            boundary edge obtained after layering.
        Much of this implementation is deferred to the following subroutines,
        which are *not* fully implemented by default:
        --> _shortenImpl()
        --> _findBoundaryMove()
        --> _setFromEdgeLocationsImpl()

        Subclasses that require this routine must therefore either:
        --> override this routine; or
        --> supply implementations for the subroutines listed above.
        In the latter case, see the documentation for these subroutines for
        details on the behaviour that must be implemented.

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
            self._setFromEdgeLocationsImpl(newEdgeLocations)

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

    def _minimiseVerticesImpl(self):
        """
        Attempts to minimise the number of vertices in the ambient
        triangulation without topologically altering this collection of
        invariant edges, potentially adding tetrahedra to do this.

        The default implementation of this routine attempts to reduce the
        number of vertices by performing some sequence of the following
        operations, subject to the constraint that this collection of
        invariant edges is not topologically altered:
        --> Shortening this collection of invariant edges using the
            _shortenImpl() routine.
        --> Close book moves, layerings and/or snap edge moves on the ambient
            triangulation.
        In particular, this implementation never creates new vertices.

        Much of this implementation is deferred to the following subroutines,
        which are *not* fully implemented by default:
        --> _shortenImpl()
        --> _minimiseBoundaryImpl()
        --> _findSnapEdge()

        Subclasses that require this routine must therefore either:
        --> override this routine; or
        --> supply implementations for the subroutines listed above.
        In the latter case, see the documentation for these subroutines for
        details on the behaviour that must be implemented.

        Adapted from Regina's Triangulation3.minimiseVertices().

        Precondition:
        --> The ambient triangulation (i.e., self.triangulation()) is valid.

        Returns:
            True if and only if a change was made to either this collection
            of invariant edges or its ambient triangulation (or both).
        """
        # Start by attempting to minimise the boundary.
        changed = self._minimiseBoundaryImpl()

        # With nothing further we can do on the boundary, all that remains is
        # to attempt to minimise the number of internal vertices.
        # We do not currently have an implementation of collapseEdge() that
        # keeps track of how edges get relabelled after the move, so we rely
        # entirely on the snapEdge() routine.
        while True:
            # Attempt shortening. This has two potential benefits:
            #   (a) Making further simplifications available later on.
            #   (b) Reducing the number of cases that we need to handle.
            if self._shortenImpl():
                changed = True

            # Is there a snap edge move available?
            snap = self._findSnapEdge()
            if snap is None:
                # Nothing further we can do.
                return changed
            changed = True
            edge, newEdgeLocations = snap

            # Perform the snap, and then update this collection of invariant
            # edges. The _findSnapEdge() routine only returns legal snap edge
            # moves, so we don't need to check again before performing.
            snapEdge( edge, False, True )
            self._setFromEdgeLocationsImpl(newEdgeLocations)
        return

    def _findSnapEdge(self):
        """
        Returns details of a snap edge move that reduces the number of
        vertices in the ambient triangulation, or None if no such moves are
        available.

        This routine only considers a snap edge move on an edge e to be
        "available" if performing this move would not topologically alter
        this collection of invariant edges. In the case where a suitable edge
        e exists, this routine returns a tuple consisting of the following:
        (0) The edge e.
        (1) A collection of edge locations that describes the location of the
            invariant edges after performing the snap edge move on the edge
            e. More precisely, after the snap edge move, the collection of
            invariant edges can be found by passing the edge locations to
            _setFromEdgeLocationsImpl().

        The InvariantEdges base class does not implement this routine, so
        subclasses that require this routine must provide an implementation.

        Pre-condition:
        --> This collection of invariant edges cannot be shortened.
        --> If the ambient triangulation has real boundary, then this
            boundary has already been minimised.
        """
        raise NotImplementedError()

    def _simplifyMonotonicImpl( self, include32 ):
        """
        Monotonically reduces the number of tetrahedra in the ambient
        triangulation without topologically altering this collection of
        invariant edges.

        The default implementation of this routine attempts to reduce the
        number of tetrahedra by performing some sequence of the following
        operations, subject to the constraint that this collection of
        invariant edges is not topologically altered:
        --> Shortening this collection of invariant edges using the
            _shortenImpl() routine.
        --> 2-0 edge moves, 2-1 edge moves and/or (optionally) 3-2 moves on
            the ambient triangulation.

        If no such simplification is possible, then this collection of
        invariant edges (and in particular, the ambient triangulation) will
        not be modified at all.

        The default implementation is partly deferred to the following
        subroutines, which are *not* implemented by default:
        --> __contains__()
        --> _updateInvariantEdges()

        Subclasses that require this routine must therefore either:
        --> override this routine; or
        --> supply implementations for the subroutines listed above.
        In the latter case, see the documentation for these subroutines for
        details on the behaviour that must be implemented.

        Each of the aforementioned elementary moves (2-0 edge, 2-1 edge and
        3-2) is performed on an edge e of the ambient triangulation, and has
        the effect of deleting e from the 1-skeleton. Thus, to ensure that
        this collection of invariant edges is not topologically altered, the
        default implementation will never perform one of these elementary
        moves on one of the invariant edges.

        The 2-0 edge and 2-1 edge moves have the additional effect of merging
        a pair of edges that share the same endpoints. If both of these
        merged edges are supposed to be invariant, then this would
        topologically alter this collection of invariant edges. However, the
        default implementation does nothing to check whether this happens;
        instead, this is deferred to the _updateInvariantEdges() subroutine.
        The subroutine is only executed *after* performing each elementary
        move, and may raise a suitable exception if it detects that the move
        topologically altered this collection of invariant edges. The default
        implementation of this routine simply passes on any such exceptions.

        Adapted from Regina's Triangulation3.simplifyToLocalMinimum() and
        SnapPea's check_for_cancellation().

        Returns:
            True if and only if the ambient triangulation was successfully
            simplified.
        """
        changed = self._shortenImpl()

        # Repeatedly attempt to reduce the size of the ambient triangulation
        # until no further reduction is possible via elementary moves.
        changedNow = True
        while True:
            changedNow = False
            for edge in self._tri.edges():
                # Make sure not to delete an invariant edge.
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
                # This move can topologically alter the collection of
                # invariant edges, but checking whether this happens is
                # deferred to the _updateInvariantEdges() routine.
                renum = twoZero(edge)
                if renum is not None:
                    changedNow = True
                    changed = True
                    break

                # Try a 2-1 edge move.
                # This move can topologically alter the collection of
                # invariant edges, but checking whether this happens is
                # deferred to the _updateInvariantEdges() routine.
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
                self._updateInvariantEdges(renum)
                self._shortenImpl()
            else:
                break

        # Nothing further we can do.
        return changed

    def _updateInvariantEdges( self, renum ):
        """
        Updates this collection of invariant edges using the given edge
        renumbering map.

        It is assumed that the renumbering map describes how the edge
        numbering of the ambient triangulation changed after performing one
        of the following elementary moves:
        --> 3-2 move
        --> 2-0 edge move
        --> 2-1 edge move
        --> 4-4 move

        The InvariantEdges base class does not implement this routine, so
        subclasses that require this routine must supply an implementation.
        Such implementations may raise an appropriate exception if the
        renumbering map comes from an elementary move that would have
        topologically altered this collection of invariant edges.
        """
        raise NotImplementedError()

    def _simplifyImpl(self):
        """
        Attempts to simplify this collection of invariant edges.

        This routine uses _minimiseVerticesImpl() and
        _simplifyMonotonicImpl(), in combination with random 4-4 moves that
        leave the invariant edges untouched.

        The default implementation is partly deferred to the following
        subroutines, which are *not* fully implemented by default:
        --> clone()
        --> _minimiseVerticesImpl()
        --> _simplifyMonotonicImpl()
        --> __contains__()
        --> _updateInvariantEdges()
        --> setAsCloneOf()

        Subclasses that require this routine must therefore either:
        --> override this routine; or
        --> supply implementations for the subroutines listed above.
        In the latter case, see the documentation for these subroutines for
        details on the behaviour that must be implemented.

        Adapted from Regina's Triangulation3.intelligentSimplify().

        Warning:
            --> Running this routine multiple times on the same loop may
                return different results, since the implementation makes
                random decisions.

        Returns:
            True if and only if this collection of invariant edges was
            successfully simplified. Otherwise, this collection of invariant
            edges will not be modified at all.
        """
        RandomEngine.reseedWithHardware()

        # Work with a clone so that we can roll back changes if we are not
        # able to reduce the number of tetrahedra.
        tempInv = self.clone()
        tempTri = tempInv.triangulation()

        # Start by minimising vertices. This will probably increase the
        # number of tetrahedra if the number of vertices is not already
        # minimal, but hopefully the monotonic simplification saves us.
        #
        # Might raise BoundsDisc.
        tempInv._minimiseVerticesImpl()
        tempInv._simplifyMonotonicImpl(True)   # Include 3-2 moves.

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
                if edge.index() in tempInv:
                    # We do not want to touch the invariant edges.
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
            fourFourChoice = fourFourAvailable[
                    RandomEngine.rand(availableCount) ]
            renum = fourFour( *fourFourChoice )
            tempInv._updateInvariantEdges(renum)
            if tempInv._simplifyMonotonicImpl(True):   # Include 3-2 moves.
                # We successfully simplified!
                # Start all over again.
                fourFourAttempts = 0
                fourFourCap = 0
            else:
                fourFourAttempts += 1

        # If simplification was successful (either by reducing the number of
        # tetrahedra, or failing that by reducing the number of vertices
        # without increasing the number of tetrahedra), then sync this
        # collection of invariant edges with the now-simpler tempInv.
        tri = self.triangulation()
        simplified = ( tempTri.size() < tri.size() or
                ( tempTri.size() == tri.size() and
                    tempTri.countVertices() < tri.countVertices() ) )
        if simplified:
            self.setAsCloneOf(tempInv)
        return simplified
