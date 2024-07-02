"""
Ideal loops for representing torus boundary components of a 3-manifold.
"""
from regina import *
from moves import threeTwo, twoZero, twoOne, fourFour


def persistentLocation(face):
    """
    Returns a persistent identifier for the location of the given face.

    In detail, this routine returns a pair (t,f), where:
    --> t is the index of a tetrahedron meeting the given face; and
    --> f is a face number of this tetrahedron that corresponds to the given
        face.
    This identifier remains valid as long as tetrahedron indices remain
    unchanged (for example, as long as no tetrahedra are ever deleted).
    """
    emb = face.embedding(0)
    return ( emb.tetrahedron().index(), emb.face() )


def snapEdge(edge):
    """
    If the endpoints of the given edge are distinct, then uses a snapped ball
    to pinch these two endpoints together.

    This operation is equivalent to performing the following two operations:
    (1) Pinching the edge, which introduces a two-tetrahedron with a single
        degree-one edge e at its heart.
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
            raise RuntimeError( "Snap edge failed unexpectedly." )
    return True


class NotLoop(Exception):
    """
    Raised when attempting to build an ideal loop from a list of edges that
    does not described an embedded closed loop.
    """
    def __init__( self, edges ):
        indices = [ e.index() for e in edges ]
        msg = ( "The edge sequence {} does not describe ".format(indices) +
                "an embedded closed loop." )
        super().__init__(msg)


class IdealLoop:
    """
    A sequence of edges representing an embedded ideal loop in a 3-manifold
    triangulation.
    """
    def __init__( self, edges=None ):
        """
        Creates an ideal loop from the given list of edges.

        If no edges are supplied, then creates an empty ideal loop with no
        data. In this case, one of the "set from" routines must be called on
        the ideal loop before performing any computations.

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

    def setFromEdges( self, edges ):
        """
        Sets this ideal loop using the given list of edges.

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

    def setFromLoop( self, loop, copyTri=True ):
        """
        Sets this ideal loop to be a clone of the given loop.

        If copyTri is True (the default), then this ideal loop will be
        embedded in a copy of the triangulation containing the given loop.
        Otherwise, if copyTri is False, then this ideal loop will be embedded
        in the same triangulation as the given loop.
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
        Sets this ideal loop using the given list of edge locations.

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
        Sets this ideal loop using a lightweight description, as constructed
        by IdealLoop.lightweightDescription().
        """
        tri = Triangulation3.fromIsoSig(sig)
        edges = []
        for tetIndex, edgeNum in edgeLocations:
            edges.append( tri.tetrahedron(tetIndex).edge(edgeNum) )
        self.setFromEdges(edges)
        return

    def _setFromRenum( self, renum ):
        """
        Sets this ideal loop using the given edge renumbering map.

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

    def triangulation(self):
        """
        Returns the triangulation that contains this ideal loop.
        """
        return self._tri

    def lightweightDescription(self):
        """
        Returns an extremely lightweight description of this ideal loop.

        In detail, this routine returns a pair (S,L), where:
        --> S is the isomorphism signature for self.triangulation(); and
        --> L is the list of edge locations in Triangulation3.fromIsoSig(S)
            corresponding to this ideal loop, with each edge location being
            given by a pair consisting of a tetrahedron index and an edge
            number.
        Thus, the returned description can be used to build an ideal loop
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

    def clone(self):
        """
        Returns a clone of this ideal loop.

        The cloned ideal loop will always be embedded in a copy of the
        triangulation containing this ideal loop.
        """
        newTri = Triangulation3( self._tri )
        newEdges = [ newTri.edge(ei) for ei in self._edgeIndices ]
        return IdealLoop(newEdges)

    def intersects( self, surf ):
        """
        Returns True if and only if this ideal loop has nonempty intersection
        with the given normal surface surf.

        Pre-condition:
        --> The given normal surface is embedded in self.triangulation().
        """
        for i in self._edgeIndices:
            if surf.edgeWeight(i).safeLongValue() > 0:
                return True
        return False

    def weight( self, surf ):
        """
        Returns the number of times this ideal loop intersects the given
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
        surface surf splits this ideal loop.

        In detail, each item of the returned list is a list of edge segments.
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

        # The given surf splits this ideal loop into multiple components, so
        # we need to do a bit more work.
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

    def minimiseVertices(self):
        """
        Reduces the number of vertices in the ambient triangulation to one.

        If the number of vertices is not already equal to one, then this
        routine increases the number of tetrahedra to achieve its goal.
        Otherwise, this routine leaves everything entirely untouched.
        """
        while self._tri.countVertices() > 1:
            # Find a suitable edge to collapse. Start with edges belonging to
            # the ideal loop, since this ensures that the loop remains
            # embedded after the collapse.
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

            # Perform the snap, and then update this ideal loop.
            snapEdge(edge)
            self.setFromEdgeLocations(edgeLocations)
        return

    def simplifyMonotonic(self):
        """
        Monotonically reduces the number of tetrahedra in the ambient
        triangulation to a local minimum, while leaving this ideal loop
        untouched.

        There should usually be no need to call this routine directly, since
        the functionality is subsumed by the more powerful simplify() and
        simplifyWithFourFour() routines.

        Adapted from Regina's Triangulation3.simplifyToLocalMinimum().

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
                # Make sure to leave the ideal loop untouched.
                if edge.index() in self._edgeIndices:
                    continue

                # Try a 3-2 move.
                renum = threeTwo(edge)
                if renum is not None:
                    changedNow = True
                    changed = True
                    break

                # Try a 2-0 edge move.
                renum = twoZero(edge)
                if renum is not None:
                    changedNow = True
                    changed = True
                    break

                # Try a 2-1 edge move.
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
                self._setFromRenum(renum)
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

        Adapted from Regina's Triangulation3.intelligentSimplify().

        Returns:
            True if and only if the ambient triangulation was successfully
            simplified. Otherwise, the ambient triangulation will not be
            modified at all.
        """
        changed = self.simplifyMonotonic()
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
            if tempLoop.simplifyMonotonic():
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

        Adapted from Regina's Triangulation3.intelligentSimplify().

        Returns:
            True if and only if this ideal loop was successfully simplified.
            Otherwise, this ideal loop will not be modified at all.
        """
        RandomEngine.reseedWithHardware()

        # Try minimising the number of vertices.
        if self._tri.countVertices() > 1:
            minVer = self.clone()
            minVer.minimiseVertices()
            minVer.simplifyWithFourFour()
            if minVer.triangulation().size() <= self._tri.size():
                self.setFromLoop( minVer, False )
                return True

        # Try minimising the number of tetrahedra.
        return self.simplifyWithFourFour()
