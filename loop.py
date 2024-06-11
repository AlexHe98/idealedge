"""
Ideal loops for representing torus boundary components of a 3-manifold.
"""
from regina import *
from moves import threeTwo, twoZero, twoOne


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
    def __init__( self, edges ):
        """
        Creates an ideal loop from the given list of edges.

        Raises NotLoop if the given list of edges does not form an embedded
        closed loop, or if the order of the edges in the given list does not
        match the order in which the edges appear in the loop.

        Pre-condition:
        --> The given list of edges is nonempty, and consists of edges that
            all belong to the same 3-manifold triangulation.
        """
        self.set(edges)

    def set( self, edges ):
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

        # Store edges of this loop in order. If we think of each edge as
        # being oriented from tail to head, so that we have an orientation on
        # the entire loop, then it's also useful to know whether vertex 0 or
        # vertex 1 is the tail.
        self._edgeIndices = []
        self._tails = []

        # Store vertices of this loop as a set, to make it as easy as
        # possible to check whether two loops are disjoint.
        self._vertIndices = set()

        # Check for degenerate case that isn't guaranteed to be ruled out by
        # the subsequent tests.
        if len(edges) == 2 and edges[0] == edges[1]:
            raise NotLoop(edges)

        # While populating the member variables, also test that the given
        # list of edges actually describes an embedded closed loop.
        firstVert = edge.vertex(0)
        lastVert = edge.vertex(0)
        broken = False
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
                broken = True
                break
        if ( broken or ( lastVert != firstVert ) or
                ( len( self._vertIndices ) != len(edges) ) ):
            raise NotLoop(edges)
        return

    def __len__(self):
        return len( self._edgeIndices )

    def __iter__(self):
        return iter( self._edgeIndices )

    def triangulation(self):
        """
        Returns the triangulation that contains this ideal loop.
        """
        return self._tri

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

    def minimiseTetrahedra(self):
        """
        Monotonically reduces the number of tetrahedra in the ambient
        triangulation to a local minimum, while leaving this ideal loop
        untouched.

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
            if not changedNow:
                break
            newEdges = []
            for ei in self._edgeIndices:
                newEdges.append( self._tri.edge( renum[ei] ) )
            self.set(newEdges)

        # Nothing further we can do.
        return changed

    def simplify(self):
        """
        Attempts to simplify this ideal loop.

        In detail, this routine attempts to reduce:
        --> the number of tetrahedra in the ambient triangulation;
        --> the number of vertices in the ambient triangulation; and
        --> the number of edges in this ideal loop.
        The priority is reducing the number of tetrahedra.

        Returns:
            True if and only if this ideal loop was successfully simplified.
            Otherwise, this ideal loop will not be modified at all.
        """
        # Work with a clone so that we can roll back changes if we fail to
        # improve the triangulation.
        use = Triangulation3( self._tri )

        # Start
        #TODO
        changedNow = True
        while changedNow:
            changedNow = False
            #TODO
            pass
        #TODO
        return
