"""
Embedded loops in a 3-manifold triangulation, which play two main roles:
--> Ideal loops for representing torus boundary components of a 3-manifold.
--> Boundary loops on triangulation with real boundary.
"""
from regina import *
from moves import twoThree, threeTwo, twoZero, twoOne, fourFour


def snapEdge(edge):
    """
    If the endpoints of the given edge are distinct, then uses a snapped ball
    to pinch these two endpoints together.

    This operation is equivalent to performing the following two operations:
    (1) Pinching the edge, which introduces a two-tetrahedron gadget with a
        single degree-one edge e at its heart.
    (2) Performing a 2-1 edge move on e.

    If the triangulation containing the given edge is currently oriented,
    then this operation will preserve the orientation.

    Pre-conditions:
    --> The given edge belongs to a triangulation with no boundary faces.

    TODO:
    --> The stated pre-conditions are stronger than necessary.

    Parameters:
    --> edge    The edge whose endpoints should be snapped together.

    Returns:
        True if and only if snapping the given edge is possible.
    """
    if edge.vertex(0) == edge.vertex(1):
        return False
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


class EmbeddedLoopException(Exception):
    pass


class NotLoop(EmbeddedLoopException):
    """
    Raised when attempting to build an embedded loop from a list of edges
    that does not describe a closed loop.
    """
    def __init__( self, edges ):
        indices = [ e.index() for e in edges ]
        msg = ( "The edge sequence {} does not describe ".format(indices) +
                "an embedded closed loop." )
        super().__init__(msg)
        return


class BoundsDisc(EmbeddedLoopException):
    """
    Raised when an embedded loop detects that it bounds an embedded disc.
    """
    def __init__(self):
        super().__init__( "Embedded loop bounds an embedded disc." )
        return


class EmbeddedLoop:
    """
    A sequence of edges representing an embedded loop in a 3-manifold
    triangulation.

    This is a base class that implements common functionality for the
    IdealLoop and BoundaryLoop classes. Although this base class can be
    instantiated, the functionality it offers is much less complete than its
    aforementioned subclasses.
    """
    def __init__( self, edges=None ):
        """
        Creates an embedded loop from the given list of edges.

        If no edges are supplied, then creates an empty object with no data.
        In this case, one of the "set from" routines must be called on the
        embedded loop before performing any computations.

        Raises NotLoop if the given list of edges does not form an embedded
        closed loop, or if the order of the edges in the given list does not
        match the order in which the edges appear in the loop.

        Pre-condition:
        --> The given list of edges is nonempty, and consists of edges that
            all belong to the same 3-manifold triangulation.
        """
        if edges is not None:
            self.setFromEdges(edges)
        return

    def setFromEdges( self, edges ):
        """
        Sets this embedded loop using the given list of edges.

        Raises NotLoop if the given list of edges does not form an embedded
        closed loop, or if the order of the edges in the given list does not
        match the order in which the edges appear in the loop.

        Pre-condition:
        --> The given list of edges is nonempty, and consists of edges that
            all belong to the same 3-manifold triangulation.
        """
        edge = edges[0]
        self._tri = edge.triangulation()
        self._reset()

        # Check for degenerate case that isn't guaranteed to be ruled out by
        # the subsequent tests.
        if len(edges) == 2 and edges[0] == edges[1]:
            raise NotLoop(edges)

        # We don't know which direction is the correct way to traverse the
        # given sequence of edges, so try both.
        firstVert = edge.vertex(0)
        lastVert = self._traverse( edges, firstVert )
        if lastVert is None:
            firstVert = edge.vertex(1)
            lastVert = self._traverse( edges, firstVert )
            if lastVert is None:
                raise NotLoop(edges)

        # Check that we actually have an embedded closed loop.
        if ( ( lastVert != firstVert ) or
                ( len( self._vertIndices ) != len(edges) ) ):
            raise NotLoop(edges)
        return

    def _reset(self):
        # Store edges of this loop in order. If we think of each edge as
        # being oriented from tail to head, so that we have an orientation on
        # the entire loop, then it's also useful to know whether vertex 0 or
        # vertex 1 is the tail.
        self._edgeIndices = []
        self._tails = []

        # Store vertices of this loop as a set, to make it as easy as
        # possible to check whether two loops are disjoint.
        self._vertIndices = set()

    def _traverse( self, edges, firstVert ):
        lastVert = firstVert
        for edge in edges:
            self._edgeIndices.append( edge.index() )
            self._vertIndices.add( lastVert.index() )

            # Find the tail of the current edge, which should join up with
            # the last vertex that we have found so far.
            if edge.vertex(0) == lastVert:
                self._tails.append(0)
                lastVert = edge.vertex(1)
            elif edge.vertex(1) == lastVert:
                self._tails.append(1)
                lastVert = edge.vertex(0)
            else:
                # Neither vertex of the current edge joins up with the last
                # vertex.
                self._reset()
                return None
        return lastVert

    def setFromLoop( self, loop, copyTri=True ):
        """
        Sets this embedded loop to be a clone of the given loop.

        If copyTri is True (the default), then this embedded loop will be
        embedded in a copy of the triangulation containing the given loop.
        Otherwise, if copyTri is False, then this embedded loop will be
        embedded in the same triangulation as the given loop.
        """
        if copyTri:
            tri = Triangulation3( loop.triangulation() )
        else:
            tri = loop.triangulation()
        edges = []
        for ei in loop:
            edges.append( tri.edge(ei) )
        self.setFromEdges(edges)
        return

    def setFromEdgeLocations( self, edgeLocations ):
        """
        Sets this embedded loop using the given list of edge locations.

        In detail, each edge location must be a pair (t, e), where t is a
        tetrahedron and e is an edge number from 0 to 5 (inclusive). Each
        tetrahedron must belong to the same triangulation.

        Raises NotLoop if the given list of edges does not form an embedded
        closed loop, or if the order of the edges in the given list does not
        match the order in which the edges appear in the loop.

        Pre-condition:
        --> The given list of edge locations is nonempty.
        --> The tetrahedra given by the first entries of each edge location
            must all belong to the same 3-manifold triangulation.
        """
        edges = []
        for tet, edgeNum in edgeLocations:
            edges.append( tet.edge(edgeNum) )
        self.setFromEdges(edges)
        return

    def setFromLightweight( self, sig, edgeLocations ):
        """
        Sets this embedded loop using a lightweight description, as
        constructed by EmbeddedLoop.lightweightDescription().
        """
        tri = Triangulation3.fromIsoSig(sig)
        edges = []
        for tetIndex, edgeNum in edgeLocations:
            edges.append( tri.tetrahedron(tetIndex).edge(edgeNum) )
        self.setFromEdges(edges)
        return

    def _setFromRenum( self, renum ):
        """
        Sets this embedded loop using the given edge renumbering map.

        This routine is for internal use only.
        """
        edges = []
        for ei in self._edgeIndices:
            edges.append( self._tri.edge( renum[ei] ) )
        self.setFromEdges(edges)
        return

    def __len__(self):
        return len( self._edgeIndices )

    def __contains__( self, edgeIndex ):
        return edgeIndex in self._edgeIndices

    def __iter__(self):
        return iter( self._edgeIndices )

    def __getitem__( self, index ):
        return self._edgeIndices[index]

    def clone(self):
        """
        Returns a clone of this embedded loop.

        The cloned loop will always be embedded in a copy of the
        triangulation containing this loop.
        """
        return EmbeddedLoop( self._cloneImpl() )

    def _cloneImpl(self):
        newTri = Triangulation3( self._tri )
        return [ newTri.edge(ei) for ei in self._edgeIndices ]

    def triangulation(self):
        """
        Returns the triangulation that contains this embedded loop.
        """
        return self._tri

    def isBoundary(self):
        """
        Does this embedded loop lie entirely in the boundary of the ambient
        triangulation?
        """
        for ei in self._edgeIndices:
            if not self._tri.edge(ei).isBoundary():
                return False
        return True

    def isInternal(self):
        """
        Does this embedded loop lie entirely in the interior of the ambient
        triangulation?
        """
        for vi in self._vertIndices:
            if self._tri.vertex(vi).isBoundary():
                return False
        return True

    def lightweightDescription(self):
        """
        Returns a lightweight description of this embedded loop.

        In detail, this routine returns a pair (S,L), where:
        --> S is the isomorphism signature for self.triangulation(); and
        --> L is the list of edge locations in Triangulation3.fromIsoSig(S)
            corresponding to this embedded loop, with each edge location
            being given by a pair consisting of a tetrahedron index and an
            edge number.
        Thus, the returned description can be used to build an embedded loop
        that is, up to combinatorial isomorphism, the same as this one.
        """
        sig, isom = self._tri.isoSigDetail()
        newEdgeLocations = []
        for ei in self._edgeIndices:
            emb = self._tri.edge(ei).embedding(0)
            oldIndex = emb.tetrahedron().index()
            newIndex = isom.simpImage(oldIndex)
            edgeNumber = Edge3.faceNumber(
                    isom.facetPerm(oldIndex) * emb.vertices() )
            newEdgeLocations.append( ( newIndex, edgeNumber ) )
        return ( sig, newEdgeLocations )

    def intersects( self, surf ):
        """
        Returns True if and only if this embedded loop has nonempty
        intersection with the given normal surface surf.

        Pre-condition:
        --> The given normal surface is embedded in self.triangulation().
        """
        for i in self._edgeIndices:
            if surf.edgeWeight(i).safeLongValue() > 0:
                return True
        return False

    def weight( self, surf ):
        """
        Returns the number of times this embedded loop intersects the given
        normal surface surf.

        Pre-condition:
        --> The given normal surface is embedded in self.triangulation().
        """
        wt = 0
        for i in self._edgeIndices:
            wt += surf.edgeWeight(i).safeLongValue()
        return wt

    def components( self, surf ):
        """
        Returns a list describing the components into which the given normal
        surface surf splits this embedded loop.

        In detail, each item of the returned list is a list of edge segments,
        where:
        --> Each list of edge segments is ordered according to the order in
            which the segments appear as we traverse this embedded loop.
        --> Each edge segment is encoded as a pair (ei, n) such that:
            --- ei is the index of the edge containing the segment in
                question; and
            --- n is the segment number (from 0 to w, inclusive, where w is
                the weight of the surface on edge ei).
        Note that for each edge e, the segments are numbered in ascending
        order from the segment incident to e.vertex(0) to the segment
        incident to e.vertex(1).
        """
        # We find all the components by simply walking around the loop. Take
        # the first component to be the one that begins *after* the first
        # point at which this loop gets split by the given surf. Thus, our
        # walk starts in the middle of the last component, so we need to make
        # sure to remember all the segments of the last component.
        lastComponent = []
        splitIndex = None
        for i in range( len(self) ):
            edgeIndex = self._edgeIndices[i]
            wt = surf.edgeWeight(edgeIndex).safeLongValue()
            if wt > 0:
                # We found the point at which the first component begins.
                if self._tails[i] == 0:
                    lastComponent.append( ( edgeIndex, 0 ) )
                    headSeg = ( edgeIndex, wt )
                else:
                    lastComponent.append( ( edgeIndex, wt ) )
                    headSeg = ( edgeIndex, 0 )
                splitIndex = i
                break
            else:
                # We are still in the middle of the last component.
                lastComponent.append( ( edgeIndex, 0 ) )
        if splitIndex is None:
            # If this loop is disjoint from the surface, then there is only
            # one component, and we have already found all the constituent
            # segments of this component.
            return [lastComponent]

        # The given surf splits this embedded loop into multiple components,
        # so we need to do a bit more work.
        components = []
        while splitIndex is not None:
            for seg in range( 1, wt ):
                # For wt >= 2, we get a sequence of short components given by
                # type-2 segments.
                components.append( [ ( edgeIndex, seg ) ] )

            # We now need to find all the segments that comprise the next
            # (long) component.
            nextComponent = [headSeg]
            continuation = splitIndex + 1

            # Unless we have already returned to the last component, we must
            # eventually find another split point at which the next component
            # begins.
            splitIndex = None
            for i in range( continuation, len(self) ):
                edgeIndex = self._edgeIndices[i]
                wt = surf.edgeWeight(edgeIndex).safeLongValue()
                if wt > 0:
                    # Found the next split point.
                    if self._tails[i] == 0:
                        nextComponent.append( ( edgeIndex, 0 ) )
                        headSeg = ( edgeIndex, wt )
                    else:
                        nextComponent.append( ( edgeIndex, wt ) )
                        headSeg = ( edgeIndex, 0 )
                    splitIndex = i
                    components.append(nextComponent)
                    break
                else:
                    # We are still in the middle of the current component.
                    nextComponent.append( ( edgeIndex, 0 ) )

        # Don't forget to include the last component.
        components.append( [ *nextComponent, *lastComponent ] )
        return components

    def _shortenImpl(self):
        # Shorten this loop by looking for triangles that intersect this loop
        # in two edges, and redirecting this loop to use the third edge.
        if len(self) < 2:
            return False
        changed = False
        redirected = True
        while redirected:
            redirected = False

            # Search for a triangle along which we can redirect.
            # Subclasses must implement the _redirectCandidates() helper to
            # find appropriate candidate triangles.
            for face in self._redirectCandidates():
                if self._attemptRedirect(face):     # Might raise BoundsDisc.
                    changed = True
                    redirected = True

                    # Back to the top.
                    break

        # No further shortening is possible.
        return changed

    def _redirectCandidates(self):
        raise NotImplementedError()

    def _attemptRedirect( self, face ):
        # If the given face intersects this loop in exactly two edges, then
        # we can redirect this loop along the third edge of the given face.
        incidentLocations = set()
        nonIncidentEdgeIndices = set()
        for e in range(3):
            edgeIndex = face.edge(e).index()
            try:
                location = self._edgeIndices.index(edgeIndex)
            except ValueError:
                # Edge is not incident to this loop.
                nonIncidentEdgeIndices.add(edgeIndex)
            else:
                # Edge is incident to this loop.
                incidentLocations.add(location)

        # Does the given face form an embedded disc bounded by this loop?
        if len(incidentLocations) == 3:
            raise BoundsDisc()

        # Perform redirect if possible.
        if len(incidentLocations) != 2:
            return False
        first, second = incidentLocations
        self._edgeIndices[first] = nonIncidentEdgeIndices.pop()
        self._edgeIndices.pop(second)
        newEdges = [ self._tri.edge(i) for i in self._edgeIndices ]
        self.setFromEdges(newEdges)
        return True

    def _minimiseBoundaryImpl(self):
        #TODO
        raise NotImplementedError()

    def _findBoundaryMove( self, bc ):
        raise NotImplementedError()


class IdealLoop(EmbeddedLoop):
    """
    A sequence of edges representing an embedded ideal loop in the interior
    of a 3-manifold triangulation.

    Some of the routines provided by this class might fail if the ideal loop
    bounds an embedded disc in the ambient triangulation (though these
    routines might nevertheless succeed in spite of the existence of such a
    disc). This class raises BoundsDisc whenever such a failure occurs.
    """
    def __init__( self, edges=None ):
        """
        Creates an ideal loop from the given list of edges.

        If no edges are supplied, then creates an empty object with no data.
        In this case, one of the "set from" routines must be called on the
        ideal loop before performing any computations.

        Raises NotLoop if the given list of edges does not form an embedded
        closed loop, or if the order of the edges in the given list does not
        match the order in which the edges appear in the loop.

        The given edges must all lie entirely in the interior of the ambient
        triangulation; in other words, after constructing the ideal loop L,
        calling L.isInternal() must return True. This condition is not
        checked, but some of the routines provided by this class might have
        undefined behaviour if this condition is not satisfied.

        Pre-condition:
        --> If supplied, the given list of edges must be nonempty, must
            consist of edges that all belong to the same 3-manifold
            triangulation T, and moreover all of these edges must lie
            entirely in the interior of T.
        """
        super().__init__(edges)
        return

    def clone(self):
        """
        Returns a clone of this ideal loop.

        The cloned ideal loop will always be embedded in a copy of the
        triangulation containing this ideal loop.
        """
        return IdealLoop( self._cloneImpl() )

    def drill(self):
        """
        Returns an ideal triangulation of the 3-manifold given by drilling
        out this loop.
        """
        drilled = Triangulation3( self._tri )
        drillLocations = []
        for ei in self._edgeIndices:
            emb = drilled.edge(ei).embedding(0)
            drillLocations.append( ( emb.tetrahedron(), emb.edge() ) )
        for tet, edgeNum in drillLocations:
            drilled.pinchEdge( tet.edge(edgeNum) )
        drilled.intelligentSimplify()
        drilled.minimiseVertices()
        drilled.intelligentSimplify()
        return drilled

    def shorten(self):
        """
        Shortens this ideal loop.

        In detail, if this ideal loop meets any triangle F in exactly two
        distinct edges, then it can be shortened by replacing these two edges
        with the third edge of F. This routine performs such shortenings
        until no further shortening is possible.

        If this ideal loop has length greater than one, then this routine
        might raise BoundsDisc.

        Returns:
            True if and only if this ideal loop was successfully shortened.
            Otherwise, this ideal loop will not be modified at all.
        """
        # Can use the base class implementation provided we supply an
        # implementation for _redirectCandidates().
        return self._shortenImpl()  # Might raise BoundsDisc.

    def _redirectCandidates(self):
        # Any triangle that is incident to this ideal loop is a candidate to
        # use for redirecting.
        for ei in self._edgeIndices:
            edge = self._tri.edge(ei)

            # Yield *all* triangles incident to current edge.
            for emb in edge.embeddings():
                yield emb.tetrahedron().triangle( emb.vertices()[3] )
            if edge.isBoundary():
                yield emb.tetrahedron().triangle( emb.vertices()[2] )
        return
    #TODO

    def minimiseBoundary(self):
        """
        Ensures that the triangulation containing this ideal loop has the
        smallest possible number of boundary triangles, potentially adding
        tetrahedra to do this.

        A side-effect of calling this routine is that it will shorten this
        ideal loop if possible.

        If this ideal loop has length greater than one, then this routine
        might raise BoundsDisc.

        The following are guaranteed to hold once this routine is finished:
        --> Every 2-sphere boundary component will have exactly two triangles
            and three vertices.
        --> Every projective plane boundary component will have exactly two
            triangles and two vertices.
        --> Every other boundary component will have exactly one vertex.

        The changes that this routine performs can always be expressed using
        only close book moves and/or layerings. In particular, this routine
        never creates new vertices, and it never creates a non-vertex-linking
        normal disc or 2-sphere if there was not one before.

        Adapted from Regina's Triangulation3.minimiseBoundary().

        Precondition:
        --> The ambient triangulation (i.e., self.triangulation()) is valid.

        Returns:
            True if the ambient triangulation was changed, or False if every
            boundary component of the ambient triangulation was already
            minimal to begin with.
        """
        # Can use the base class implementation provided we supply an
        # implementation for _findBoundaryMove().
        return self._minimiseBoundaryImpl()

    def _findBoundaryMove( self, bc ):
        # Precondition:
        #   --> This loop cannot be shortened
        #   --> bc is a boundary component of self._tri, and is not minimal

        # First try to find a close book move, which does not increase the
        # number of tetrahedra.
        #TODO

        # We could not find a close book move.
        # In this case, because bc is supposed to be non-minimal, there must
        # be a boundary edge e that joins two distinct vertices, and we can
        # simplify bc by layering across e and then performing a close book
        # move on the newly layered edge. (This is equivalent to attaching a
        # snapped ball along e.)
        #TODO
        raise NotImplementedError()

    def minimiseVertices(self):
        """
        Reduces the number of vertices in the ambient triangulation to one.

        If the number of vertices is not already equal to one, then this
        routine increases the number of tetrahedra to achieve its goal.
        Otherwise, this routine leaves everything entirely untouched.

        If this ideal loop has length greater than one, then this routine
        might raise BoundsDisc.
        """
        while self._tri.countVertices() > 1:
            # Might raise BoundsDisc.
            self.shorten()

            # Find a suitable edge to collapse. Start with edges belonging to
            # the ideal loop.
            if len(self) > 1:
                edge = self._tri.edge( self._edgeIndices[0] )
            else:
                for edge in self._tri.edges():
                    if edge.vertex(0) != edge.vertex(1):
                        break

            # Make sure we'll be able to find the new ideal loop after
            # snapping the endpoints of this edge together.
            edgeLocations = []
            for ei in self._edgeIndices:
                if ei == edge.index():
                    continue
                emb = self._tri.edge(ei).embedding(0)
                edgeLocations.append( ( emb.tetrahedron(), emb.edge() ) )

            # Perform the snap, and then update this ideal loop. Note that
            # this will not have any unintended side-effects on the ideal
            # loop, both because we ran shorten() and because we chose to
            # prioritise collapsing edges belonging to the ideal loop.
            snapEdge(edge)
            self.setFromEdgeLocations(edgeLocations)
        return

    def simplifyBasic(self):
        """
        Uses 2-0 and 2-1 edge moves to reduce the number of tetrahedra in the
        ambient triangulation, while leaving this ideal loop untouched.

        There should usually be no need to call this routine directly, since
        the functionality is subsumed by the more powerful simplify(),
        simplifyWithFourFour() and simplifyMonotonic() routines.

        If this ideal loop has length greater than one, then this routine
        might raise BoundsDisc.

        Adapted from SnapPea's check_for_cancellation().

        Returns:
            True if and only if the ambient triangulation was successfully
            simplified. Otherwise, the ambient triangulation will not be
            modified at all.
        """
        # Do not include 3-2 moves.
        # Might raise BoundsDisc.
        return self._simplifyImpl(False)

    def simplifyMonotonic(self):
        """
        Monotonically reduces the number of tetrahedra in the ambient
        triangulation to a local minimum, while leaving this ideal loop
        untouched.

        There should usually be no need to call this routine directly, since
        the functionality is subsumed by the more powerful simplify() and
        simplifyWithFourFour() routines.

        If this ideal loop has length greater than one, then this routine
        might raise BoundsDisc.

        Adapted from Regina's Triangulation3.simplifyToLocalMinimum().

        Returns:
            True if and only if the ambient triangulation was successfully
            simplified. Otherwise, the ambient triangulation will not be
            modified at all.
        """
        # Include 3-2 moves.
        # Might raise BoundsDisc.
        return self._simplifyImpl(True)

    def _simplifyImpl( self, include32 ):
        changed = False     # Has anything changed ever?    (Return value.)
        changedNow = True   # Did we just change something? (Loop control.)
        while True:
            changedNow = False
            for edge in self._tri.edges():
                # Make sure to leave the ideal loop untouched.
                if edge.index() in self._edgeIndices:
                    continue

                # If requested, try a 3-2 move.
                if include32:
                    renum = threeTwo(edge)
                    if renum is not None:
                        changedNow = True
                        changed = True
                        break

                # Try a 2-0 edge move.
                # This move can destroy the ideal loop if it bounds a disc.
                renum = twoZero(edge)
                if renum is not None:
                    changedNow = True
                    changed = True
                    break

                # Try a 2-1 edge move.
                # This move can destroy the ideal loop if it bounds a disc.
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
            # the details of the ideal loop, and then check whether we can
            # improve further.
            if changedNow:
                try:
                    # If we destroyed the ideal loop, then this will raise
                    # NotLoop.
                    self._setFromRenum(renum)
                except NotLoop:
                    # As noted above, the ideal loop can only get destroyed
                    # if it bounds a disc.
                    raise BoundsDisc()
            else:
                break

        # Nothing further we can do.
        return changed

    def simplifyWithFourFour(self):
        """
        With the assistance of random 4-4 moves, attempts to reduce the
        number of tetrahedra in the ambient triangulation, while leaving this
        ideal loop untouched.

        In more detail, this routine uses simplifyMonotonic() until it
        reaches a local minimum, and then uses random 4-4 moves to attempt to
        escape this local minimum.

        Unlike the simplify() routine, this routine never changes the number
        of vertices.

        If this ideal loop has length greater than one, then this routine
        might raise BoundsDisc.

        Adapted from Regina's Triangulation3.intelligentSimplify().

        Returns:
            True if and only if the ambient triangulation was successfully
            simplified. Otherwise, the ambient triangulation will not be
            modified at all.
        """
        changed = self.simplifyMonotonic()  # Might raise BoundsDisc.
        tempLoop = self.clone()
        tempTri = tempLoop.triangulation()

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
                    # We do not want to touch the ideal loop.
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
            tempLoop._setFromRenum(renum)
            if tempLoop.simplifyMonotonic():    # Might raise BoundsDisc.
                # We successfully simplified!
                # Start all over again.
                fourFourAttempts = 0
                fourFourCap = 0
            else:
                fourFourAttempts += 1

        # If the 4-4 moves were successful, then sync this ideal loop with
        # the now-simplified tempLoop.
        if tempTri.size() < self.triangulation().size():
            self.setFromLoop(tempLoop)
            changed = True
        return changed

    def simplify(self):
        """
        Attempts to simplify this ideal loop.

        The rationale of this routine is to combine the functionality of the
        minimiseVertices() and simplifyWithFourFour() routines. In detail,
        this routine attempts to reduce the number of vertices in the ambient
        triangulation to one without increasing the number of tetrahedra. If
        this fails, or if the number of vertices is already equal to one,
        then this routine attempts to minimise the number of tetrahedra
        without altering the number of vertices.

        If this ideal loop has length greater than one, then this routine
        might raise BoundsDisc.

        Adapted from Regina's Triangulation3.intelligentSimplify().

        Returns:
            True if and only if this ideal loop was successfully simplified.
            Otherwise, this ideal loop will not be modified at all.
        """
        RandomEngine.reseedWithHardware()

        # Try minimising the number of vertices.
        if self._tri.countVertices() > 1:
            minVer = self.clone()
            minVer.minimiseVertices()       # Might raise BoundsDisc.
            minVer.simplifyWithFourFour()   # Might raise BoundsDisc.
            if minVer.triangulation().size() <= self._tri.size():
                self.setFromLoop( minVer, False )
                return True

        # Try minimising the number of tetrahedra.
        return self.simplifyWithFourFour()  # Might raise BoundsDisc.

    def randomise(self):
        """
        Attempts to randomly retriangulate this ideal loop.

        If this ideal loop has length greater than one, then this routine
        might raise BoundsDisc.

        Adapted from SnapPea's randomize_triangulation().
        """
        RandomEngine.reseedWithHardware()
        randomisation = 4       # Hard-coded value copied from SnapPea.
        count = randomisation * self._tri.size()
        while count > 0:
            count -= 1

            # Attempt a random 2-3 move.
            renum = twoThree( self._tri.triangle(
                RandomEngine.rand( self._tri.countTriangles() ) ) )
            if renum is not None:
                self._setFromRenum(renum)

                # Try to force future random 2-3 moves to make "interesting"
                # changes.
                self.simplifyBasic()    # Might raise BoundsDisc.

        # Finish up by simplifying. The built-in randomness should hopefully
        # take us somewhere new.
        self.simplify()     # Might raise BoundsDisc.
        return


class BoundaryLoop(EmbeddedLoop):
    """
    A sequence of edges representing an embedded loop on the boundary of a
    3-manifold triangulation.

    Some of the routines provided by this class might fail if the boundary
    loop bounds an embedded disc in the ambient triangulation (though these
    routines might nevertheless succeed in spite of the existence of such a
    disc). This class raises BoundsDisc whenever such a failure occurs.
    """
    def __init__( self, edges=None ):
        """
        Creates a boundary loop from the given list of edges.

        If no edges are supplied, then creates an empty object with no data.
        In this case, one of the "set from" routines must be called on the
        boundary loop before performing any computations.

        Raises NotLoop if the given list of edges does not form an embedded
        closed loop, or if the order of the edges in the given list does not
        match the order in which the edges appear in the loop.

        The given edges must all lie entirely on the boundary of the ambient
        triangulation; in other words, after constructing the boundary loop
        L, calling L.isBoundary() must return True. This condition is not
        checked, but some of the routines provided by this class might have
        undefined behaviour if this condition is not satisfied.

        Pre-condition:
        --> If supplied, the given list of edges must be nonempty, must
            consist of edges that all belong to the same 3-manifold
            triangulation T, and moreover all of these edges must lie
            entirely on the boundary of T.
        """
        super().__init__(edges)
        return

    def clone(self):
        """
        Returns a clone of this boundary loop.

        The cloned boundary loop will always be embedded in a copy of the
        triangulation containing this boundary loop.
        """
        return BoundaryLoop( self._cloneImpl() )

    def shorten(self):
        """
        Shortens this boundary loop.

        In detail, if this boundary loop meets any *boundary* triangle F in
        exactly two distinct edges, then it can be shortened by replacing
        these two edges with the third edge of F. This routine performs such
        shortenings until no further shortening is possible.

        If this boundary loop has length greater than one, then this routine
        might raise BoundsDisc.

        Returns:
            True if and only if this boundary loop was successfully
            shortened. Otherwise, this boundary loop will not be modified at
            all.
        """
        # Can use the base class implementation provided we supply an
        # implementation for _redirectCandidates().
        return self._shortenImpl()  # Might raise BoundsDisc.

    def _redirectCandidates(self):
        # For boundary loops, we only want to redirect along *boundary*
        # triangles (to ensure that the loop stays in the boundary).
        for ei in self._edgeIndices:
            edge = self._tri.edge(ei)

            # Yield only *boundary* triangles incident to current edge.
            # (Note that as a precondition, the current edge is assumed to be
            # a boundary edge.)
            front = edge.front()
            yield front.tetrahedron().triangle( front.vertices()[3] )
            back = edge.back()
            yield back.tetrahedron().triangle( back.vertices()[2] )
        return
    #TODO

    def minimiseBoundary(self):
        """
        Ensures that the triangulation containing this boundary loop has the
        smallest possible number of boundary triangles, potentially adding
        tetrahedra to do this.

        If this boundary loop has length greater than one, then this routine
        might raise BoundsDisc.

        The following are guaranteed to hold once this routine is finished:
        --> Every 2-sphere boundary component will have exactly two triangles
            and three vertices.
        --> Every projective plane boundary component will have exactly two
            triangles and two vertices.
        --> Every other boundary component will have exactly one vertex.

        The changes that this routine performs can always be expressed using
        only close book moves and/or layerings. In particular, this routine
        never creates new vertices, and it never creates a non-vertex-linking
        normal disc or 2-sphere if there was not one before.

        Adapted from Regina's Triangulation3.minimiseBoundary().

        Precondition:
        --> The ambient triangulation (i.e., self.triangulation()) is valid.

        Returns:
            True if the ambient triangulation was changed, or False if every
            boundary component of the ambient triangulation was already
            minimal to begin with.
        """
        # Can use the base class implementation provided we supply an
        # implementation for _findBoundaryMove().
        return self._minimiseBoundaryImpl()

    def _findBoundaryMove( self, bc ):
        # Precondition:
        #   --> This loop cannot be shortened
        #   --> bc is a boundary component of self._tri, and is not minimal

        # Prioritise reducing the length of this loop, and try to do so
        # without adding too many tetrahedra.
        #TODO
        raise NotImplementedError()
