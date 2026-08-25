"""
Embedded loops in a 3-manifold triangulation, which play two main roles:
--> Ideal loops for representing torus boundary components of a 3-manifold.
--> Boundary loops on triangulations with real boundary.
"""
from regina import *
from aux.looperror import NotLoop, BoundsDisc
from aux.edgeemb import edgesFromEmbeddings, edgeOrientationFromEmbedding
from segment import OrientedSegment
from chord import NormalChord


class EmbeddedLoop:
    """
    A sequence of edges representing an embedded loop in a 3-manifold
    triangulation.

    This is a base class that implements common functionality for the
    IdealLoop and BoundaryLoop classes. Although this base class can be
    instantiated, the functionality it offers is much less complete than its
    aforementioned subclasses.

    A core feature of this class is that it effectively stores a list of edge
    *indices* corresponding to the edges of the embedded loop; moreover, the
    order that these edge indices appear in the list corresponds to an
    orientation on the loop. Thus, for any instance loop of this class, the
    following functionality is available:
    --> (e in loop) is True if and only if loop.triangulation().edge(e) is an
        edge in the loop
    --> len(loop) is the number of edges in the loop
    --> iterating through the loop yields all the edge indices in an order
        that matches the loop's orientation
    --> for i between 0 and (len(loop) - 1), inclusive, loop[i] returns the
        index of the ith edge in the loop, and again the order matches the
        loop's orientation
    """
    def __init__( self, edges=None, orientation=0 ):
        """
        Creates an embedded loop from the given list of edges.

        If no edges are supplied, then creates an empty object with no data.
        In this case, one of the "set from" routines must be called on the
        embedded loop before performing any computations.

        If the optional orientation argument is not supplied, then the
        embedded loop will be assigned an arbitrary orientation. Otherwise,
        the supplied orientation must be one of the following:
        --> +1 if the first edge in the given list of edges should be oriented
            oriented from vertex 0 to vertex 1;
        --> -1 if the first edge should be oriented from vertex 1 to vertex 0;
            or
        --> 0 if this routine should be allowed to choose an arbitrary
            orientation.

        Raises NotLoop if the given list of edges does not form an embedded
        closed loop, or if the order of the edges in the given list does not
        match the order in which the edges appear in the loop.

        Precondition:
        --> The given list of edges is nonempty, and consists of edges that
            all belong to the same 3-manifold triangulation.
        """
        #NOTE For reference, it is often helpful to work with orientations
        #   using edge embeddings. In detail, let e = edges[0], and let
        #   emb = e.embedding(i) for any i in { 0, ..., e.degree() - 1 }
        #   (everything we say here will be independent of the choice of i).
        #
        #   Orientation +1 means that e is embedded in emb.tetrahedron() with
        #   an orientation that runs away from emb.vertices()[0] and towards
        #   emb.vertices()[1]. On the other hand, orientation -1 means that
        #   the orientation on e runs away from emb.vertices()[1] and towards
        #   emb.vertices()[0]. Pictorially, the diagram below left shows a +1
        #   orientation, and the diagram below right shows a -1 orientation.
        #
        #         emb.vertices()[1]          emb.vertices()[1]
        #                 •                          •
        #                /|\                        /|\
        #               / | \                      / | \
        #              /  ↑  \                    /  ↓  \
        #             •---|---•                  •---|---•
        #              \  ↑  /                    \  ↓  /
        #               \ | /                      \ | /
        #                \|/                        \|/
        #                 •                          •
        #         emb.vertices()[0]          emb.vertices()[0]
        #
        if edges is not None:
            self.setFromEdges( edges, orientation )
        return

    def setFromEdges( self, edges, orientation=0 ):
        """
        Sets this embedded loop using the given list of edges.

        If the optional orientation argument is not supplied, then the
        embedded loop will be assigned an arbitrary orientation. Otherwise,
        the supplied orientation must be one of the following:
        --> +1 if the first edge in the given list of edges should be oriented
            oriented from vertex 0 to vertex 1;
        --> -1 if the first edge should be oriented from vertex 1 to vertex 0;
            or
        --> 0 if this routine should be allowed to choose an arbitrary
            orientation.

        Raises NotLoop if the given list of edges does not form an embedded
        closed loop, or if the order of the edges in the given list does not
        match the order in which the edges appear in the loop.

        Precondition:
        --> The given list of edges is nonempty, and consists of edges that
            all belong to the same 3-manifold triangulation.
        """
        if orientation not in {1,-1,0}:
            raise ValueError(
                    "If supplied, orientation must be either +1, -1 or 0." )
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

        # If necessary, reverse orientation to match the user specification.
        if ( ( orientation is not None ) and
            ( self.orientation() != orientation ) ):
            self.reverseOrientation()
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

    def __str__(self):
        return "Edge indices {}, Orientation {}".format(
                self._edgeIndices, self.orientation() )

    def __repr__(self):
        return "<loop.{}: {}>".format( type(self).__name__, str(self) )

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

    @classmethod
    def fromEdgeEmbeddings( cls, edgeEmbeddings, orientation=0 ):
        """
        Constructs an embedded loop using the given list of edge embeddings.

        The elements in edgeEmbeddings must all be EdgeEmbedding3 objects
        referencing tetrahedra in the same Triangulation3 object.

        If the optional orientation argument is not supplied, then the
        embedded loop will be assigned an arbitrary orientation. Otherwise,
        the supplied orientation must be one of the following:
        --> +1 if the first edge described by the given list of edge
            embeddings should be oriented from vertex 0 to vertex 1 (here,
            vertex numbers are with respect to the edge embedding, which might
            differ from the vertex numbers of the underlying edge if the
            ambient triangulation has been modified since the edge embedding
            was constructed);
        --> -1 if the first edge should be oriented from vertex 1 to vertex 0;
            or
        --> 0 if this routine should be allowed to choose an arbitrary
            orientation.

        Raises NotLoop if the given list of edges does not form an embedded
        closed loop, or if the order of the edges in the given list does not
        match the order in which the edges appear in the loop.

        Precondition:
        --> The given list of edge embeddings is nonempty.
        --> The given edge embeddings must all reference tetrahedra belonging
            to the same 3-manifold triangulation.
        """
        return cls( edgesFromEmbeddings(edgeEmbeddings),
                   edgeOrientationFromEmbedding(
                       edgeEmbeddings[0], orientation ) )

    def _setFromRelab( self, relab ):
        """
        Sets this embedded loop using the relabelling described by the given
        EdgeLabelling.

        This routine is for internal use only. The purpose of this routine is
        to update the embedded loop whenever the ambient triangulation has
        been modified by a local move. See the twoThree, threeTwo, twoZero,
        twoOne, and fourFour routines from retriangulate/moves.py for examples
        of how relabellings are specified.

        Pre-condition:
        --> The given EdgeLabelling relab tracks every index ei in self.
        """
        edges = []
        oldOrientation = self.edgeOrientation(0)
        newOrientation = 0
        for ei in self:
            edge = self._tri.edge(
                    relab.underlyingEdgeIndex(ei) )
            edges.append(edge)
            if newOrientation == 0:
                # We are looking at edge 0 of the loop, which is the edge that
                # determines the orientation of the loop. The embedding
                # relab[ei] will orient this edge in the same direction as
                # before, whereas the underlying labelling of the edge might
                # or might not be the same as before. We therefore need to
                # compare relab[ei] with the underlying labelling to determine
                # the newOrientation for the loop.
                referenceEmb = relab[ei]
                vertexPerm = referenceEmb.tetrahedron().edgeMapping(
                        referenceEmb.edge() )
                if vertexPerm[0] == referenceEmb.vertices()[0]:
                    newOrientation = oldOrientation
                else:
                    newOrientation = -1 * oldOrientation
        self.setFromEdges( edges, newOrientation )
        return

    def __len__(self):
        return len( self._edgeIndices )

    def __contains__( self, edgeIndex ):
        return edgeIndex in self._edgeIndices

    def __iter__(self):
        return iter( self._edgeIndices )

    def __getitem__( self, index ):
        return self._edgeIndices[index]

    def vertexIndices(self):
        """
        Yields all indices of vertices incident to this embedded loop (in no
        particular order).
        """
        return iter( self._vertIndices )

    def edgeIndices(self):
        """
        Returns (a copy of) the underlying list of edge indices in this
        embedded loop.
        """
        return list( self._edgeIndices )

    def clone(self):
        """
        Returns a clone of this embedded loop.

        The cloned loop will always be embedded in a copy of the
        triangulation containing this loop.
        """
        # We use the built-in type() function to make sure that subclasses
        # will construct clones of the correct type.
        return type(self)( self._cloneImpl(), self.orientation() )

    def _cloneImpl( self, newTri=None ):
        """
        Returns a list of edges which can be used to create a clone of this
        EmbeddedLoop.

        If newTri is None (the default), then the edges in the returned list
        will all belong to a newly-constructed copy of self.triangulation().
        Otherwise, newTri must be a triangulation that is combinatorially
        identical to self.triangulation().
        """
        if newTri is None:
            newTri = Triangulation3( self._tri )
        return [ newTri.edge(ei) for ei in self ]

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
        for ei in self:
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

    def intersects( self, surf ):
        """
        Returns True if and only if this embedded loop has nonempty
        intersection with the given normal surface surf.

        Precondition:
        --> The given normal surface is embedded in self.triangulation().
        """
        for i in self:
            if surf.edgeWeight(i).pythonValue() > 0:
                return True
        return False

    def weight( self, surf ):
        """
        Returns the number of times this embedded loop intersects the given
        normal surface surf.

        Precondition:
        --> The given normal surface is embedded in self.triangulation().
        """
        wt = 0
        for i in self:
            wt += surf.edgeWeight(i).pythonValue()
        return wt

    #TODO Want a way to test disjointness of pairs of loops. A generalisable
    #   solution would be a commonVertices(otherLoop) routine, but do we want
    #   something whose name more obviously corresponds to this use case?

    def splitIntoChords( self, surf ):
        """
        Splits this embedded loop along the given normal surface, and returns
        a set containing the resulting NormalChord objects.

        Each NormalChord in the returned set will be oriented in the same
        direction as this loop.

        If this embedded loop is disjoint from surf, then the returned set
        will of course contain just a single chord. The two ends of this
        chord will be abstractly joined with each other to indicate that no
        actual split occurred along surf.
        """
        # We find all the chords by simply walking around the loop. Take the
        # first chord to be the one that begins *after* the first point at
        # which this loop gets split by the given surf. Thus, our walk starts
        # in the middle of the last chord, so we need to make sure to
        # remember all the segments of the last chord.
        lastChordSegs = []
        splitIndex = None
        for i in range( len(self) ):
            edgeIndex = self._edgeIndices[i]
            orientation = self.edgeOrientation(i)
            wt = surf.edgeWeight(edgeIndex).pythonValue()
            tailSeg = OrientedSegment(
                    surf, edgeIndex, 0, orientation )
            if wt > 0:
                # We found the point at which the first chord begins.
                headSeg = OrientedSegment(
                        surf, edgeIndex, wt, orientation )
                if self._tails[i] == 1:
                    tailSeg, headSeg = headSeg, tailSeg
                lastChordSegs.append(tailSeg)
                splitIndex = i
                break
            else:
                # We are still in the middle of the last chord.
                lastChordSegs.append(tailSeg)

        if splitIndex is None:
            # If this loop is disjoint from the surface, then there is only
            # one chord, and we have already found all the constituent
            # segments of this chord. All that remains is to abstractly join
            # the two ends back together.
            onlyChord = NormalChord(lastChordSegs)
            onlyChord.join( 0, onlyChord, 1 )
            return {onlyChord}

        # The given surf splits this embedded loop into multiple chords, so
        # we need to do a bit more work.
        chordSet = set()
        while splitIndex is not None:
            # Abuse the fact that the following variables all persist beyond
            # the scope of the above for loop.
            #   --> headSeg
            #   --> wt
            #   --> edgeIndex
            orientation = headSeg.orientation()
            for segPos in range( 1, wt ):
                # For wt >= 2, we get a sequence of short chords given by
                # type-2 segments.
                chordSet.add( NormalChord( [ OrientedSegment(
                    surf, edgeIndex, segPos, orientation ) ] ) )

            # We now need to find all the segments that comprise the next
            # (long) chord.
            nextChordSegs = [headSeg]
            continuation = splitIndex + 1

            # Unless we have already returned to the last chord, we must
            # eventually find another split point at which the next chord
            # begins.
            splitIndex = None
            for i in range( continuation, len(self) ):
                edgeIndex = self._edgeIndices[i]
                orientation = self.edgeOrientation(i)
                wt = surf.edgeWeight(edgeIndex).pythonValue()
                tailSeg = OrientedSegment(
                        surf, edgeIndex, 0, orientation )
                if wt > 0:
                    # Found the next split point.
                    headSeg = OrientedSegment(
                            surf, edgeIndex, wt, orientation )
                    if self._tails[i] == 1:
                        tailSeg, headSeg = headSeg, tailSeg
                    nextChordSegs.append(tailSeg)
                    splitIndex = i
                    chordSet.add( NormalChord(nextChordSegs) )
                    break
                else:
                    # We are still in the middle of the current chord.
                    nextChordSegs.append(tailSeg)

        # Don't forget to include the last chord.
        chordSet.add( NormalChord( [ *nextChordSegs, *lastChordSegs ] ) )
        assert ( sum( [ len(chord) for chord in chordSet ] ) ==
                len(self) + self.weight(surf) )
        return chordSet

    def orientation(self):
        """
        Returns the orientation of this embedded loop.

        In detail, the returned value will be:
        --> +1 if the first edge e in this embedded loop is oriented from
            vertex 0 to vertex 1; and
        --> -1 if e is oriented from vertex 1 to vertex 0.
        """
        return self.edgeOrientation(0)

    def edgeOrientation( self, index ):
        """
        Returns the orientation of the edge at the given index in this
        embedded loop.

        In detail, the returned value will be:
        --> +1 if the edge is oriented from vertex 0 to vertex 1; and
        --> -1 if the edge is oriented from vertex 1 to vertex 0.
        """
        if self._tails[index] == 0:
            # Tail at 0, head at 1.
            return 1
        # Tail at 1, head at 0.
        return -1

    def reverseOrientation(self):
        """
        Reverses the orientation of this embedded loop.

        This routine preserves the first edge in this embedded loop, and
        therefore ensures that the sign of self.orientation() is actually
        reversed after calling this routine.
        """
        # Preserve the first edge, but reverse the order in which we traverse
        # the other edges of this loop.
        self._edgeIndices[1:] = self._edgeIndices[:0:-1]
        self._tails[1:] = self._tails[:0:-1]

        # Switch head and tail for each edge of this loop.
        self._tails = [ 1 - i for i in self._tails ]
        return

    def shorten(self):
        """
        Shortens this embedded loop by looking for triangles that intersect
        this loop in two edges, and redirecting this loop to use the third
        edge.

        The default implementation of this routine requires the helper
        routine _redirectCandidates(), which is *not* implemented by default.
        Thus, subclasses that require this routine must either:
        --> override this routine; or
        --> supply an implementation for _redirectCandidates().
        In the latter case, see the documentation for _redirectCandidates()
        for details on the behaviour that must be implemented.

        If this loop bounds a disc, then this routine might (but is not
        guaranteed to) raise BoundsDisc.

        The default implementation raises BoundsDisc if and only if the
        _redirectCandidates() routine yields a face that is incident to this
        loop in three edges (in such a case, the face forms an embedded disc
        with boundary given by this loop).

        Returns:
            True if and only if this embedded loop was successfully
            shortened. In the case where no shortening occurred, this
            embedded loop will remain entirely untouched.
        """
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
        """
        Yields candidate triangles of self.triangulation() across which the
        shorten() routine should attempt to redirect this loop.

        The EmbeddedLoop base class does not implement this routine, so
        subclasses that require this routine must provide an implementation.
        """
        raise NotImplementedError()

    def _attemptRedirect( self, face ):
        r"""
        Attempts to redirect this loop across the given face.

        If the given face intersects this loop in exactly two *distinct*
        edges, then we can redirect this loop along the third edge of the
        given face.

                    Before redirect         After redirect
                           •                       •
                          / \
                         /   \
                        /     \
                       •       •               •-------•

        If such a redirect is possible, then this routine performs the
        redirect and returns True. Otherwise, this routine leaves this loop
        entirely untouched and returns False.

        If this loop is incident to all three edges of the given face, then
        the face forms an embedded disc with boundary given by this loop. In
        such a situation, this routine raises BoundsDisc.

        Parameters
        --> face    the triangular face across which to attempt a redirect of
                    this loop

        Returns:
            True if and only if this routine successfully performs the
            requested redirect
        """
        incidentLocations = set()       # Location in this loop.
        tails = set()
        heads = set()
        nonIncidentEdgeNums = set()
        for e in range(3):
            try:
                location = self._edgeIndices.index(
                        face.edge(e).index() )
            except ValueError:
                nonIncidentEdgeNums.add(e)
            else:
                incidentLocations.add(location)

                # Which vertices of this face are at the tail and head?
                eTail = self._tails[location]
                tails.add( face.edgeMapping(e)[eTail] )
                heads.add( face.edgeMapping(e)[1-eTail] )

        # Does the given face form an embedded disc bounded by this loop?
        # Equivalently, is the given face incident to this loop in 3 distinct
        # edges?
        if len(incidentLocations) == 3:
            raise BoundsDisc()

        # Perform redirect if possible.
        if len(incidentLocations) != 2:
            return False
        swap, delete = incidentLocations
        newEdgeNum = nonIncidentEdgeNums.pop()
        if ( swap == 0 ) or ( delete == 0 and swap == 1 ):
            # New orientation is determined by the new edge.
            #
            # After redirecting, the tail of the new edge is given by the old
            # tail that is not also an old head.
            newTail = face.edgeMapping(newEdgeNum).inverse()[
                    (tails - heads).pop() ]
            if newTail == 0:
                newOrientation = 1
            else:
                newOrientation = -1
        else:
            # New orientation is determined by an edge that is already part
            # of this loop.
            if delete == 0:
                newOrientation = self.edgeOrientation(1)
            else:
                newOrientation = self.edgeOrientation(0)
        self._edgeIndices[swap] = face.edge(newEdgeNum).index()
        self._edgeIndices.pop(delete)
        newEdges = [ self._tri.edge(i) for i in self ]
        self.setFromEdges( newEdges, newOrientation )
        return True


class IdealLoop(EmbeddedLoop):
    """
    A sequence of edges representing an embedded ideal loop in the interior
    of a 3-manifold triangulation.

    Some of the routines provided by this class might fail if the ideal loop
    bounds an embedded disc in the ambient triangulation (though these
    routines might nevertheless succeed in spite of the existence of such a
    disc). This class raises BoundsDisc whenever such a failure occurs.

    A core feature of this class is that it effectively stores a list of edge
    *indices* corresponding to the edges of the ideal loop; moreover, the
    order that these edge indices appear in the list corresponds to an
    orientation on the loop. Thus, for any instance loop of this class, the
    following functionality is available:
    --> (e in loop) is True if and only if loop.triangulation().edge(e) is an
        edge in the loop
    --> len(loop) is the number of edges in the loop
    --> iterating through the loop yields all the edge indices in an order
        that matches the loop's orientation
    --> for i between 0 and (len(loop) - 1), inclusive, loop[i] returns the
        index of the ith edge in the loop, and again the order matches the
        loop's orientation
    """
    def __init__( self, edges=None, orientation=0 ):
        """
        Creates an ideal loop from the given list of edges.

        If no edges are supplied, then creates an empty object with no data.
        In this case, one of the "set from" routines must be called on the
        ideal loop before performing any computations.

        If the optional orientation argument is not supplied, then the
        embedded loop will be assigned an arbitrary orientation. Otherwise,
        the supplied orientation must be one of the following:
        --> +1 if the first edge in the given list should be oriented from
            vertex 0 to vertex 1;
        --> -1 if the first edge should be oriented from vertex 1 to vertex 0;
            or
        --> 0 if this routine should be allowed to choose an arbitrary
            orientation.

        Raises NotLoop if the given list of edges does not form an embedded
        closed loop, or if the order of the edges in the given list does not
        match the order in which the edges appear in the loop.

        The given edges must all lie entirely in the interior of the ambient
        triangulation; in other words, after constructing the ideal loop L,
        calling L.isInternal() must return True. This condition is not
        checked, but some of the routines provided by this class might have
        undefined behaviour if this condition is not satisfied.

        Precondition:
        --> If supplied, the given list of edges must be nonempty, must
            consist of edges that all belong to the same 3-manifold
            triangulation T, and moreover all of these edges must lie
            entirely in the interior of T.
        """
        super().__init__( edges, orientation )
        return

    def shorten(self):
        """
        Shortens this ideal loop.

        In detail, if this ideal loop meets any triangle F in exactly two
        distinct edges, then it can be shortened by replacing these two edges
        with the third edge of F. This routine performs such shortenings
        until no further shortening is possible.

        This routine raises BoundsDisc if self.triangulation() includes a
        triangular face F that forms an embedded disc whose boundary is given
        by this ideal loop.

        Returns:
            True if and only if this ideal loop was successfully shortened.
            In the case where no shortening occurred, this ideal loop will
            remain entirely untouched.
        """
        # Can use the default implementation provided we supply an
        # implementation for _redirectCandidates().
        # Since _redirectCandidates() yields all triangles incident to this
        # loop, this will raise BoundsDisc as promised in the documentation.
        return super().shorten()

    def _redirectCandidates(self):
        """
        Yields candidate triangles of self.triangulation() across which the
        shorten() routine should attempt to redirect this loop.

        For an IdealLoop, every triangle incident to the loop is a candidate.
        """
        for ei in self:
            edge = self._tri.edge(ei)

            # Yield *all* triangles incident to current edge.
            # Note that as a precondition, the edge is assumed to be internal.
            for emb in edge.embeddings():
                yield emb.tetrahedron().triangle( emb.vertices()[3] )
        return


class BoundaryLoop(EmbeddedLoop):
    """
    A sequence of edges representing an embedded loop on the boundary of a
    3-manifold triangulation.

    Some of the routines provided by this class might fail if the boundary
    loop bounds an embedded disc in the ambient triangulation (though these
    routines might nevertheless succeed in spite of the existence of such a
    disc). This class raises BoundsDisc whenever such a failure occurs.

    A core feature of this class is that it effectively stores a list of edge
    *indices* corresponding to the edges of the boundary loop; moreover, the
    order that these edge indices appear in the list corresponds to an
    orientation on the loop. Thus, for any instance loop of this class, the
    following functionality is available:
    --> (e in loop) is True if and only if loop.triangulation().edge(e) is an
        edge in the loop
    --> len(loop) is the number of edges in the loop
    --> iterating through the loop yields all the edge indices in an order
        that matches the loop's orientation
    --> for i between 0 and (len(loop) - 1), inclusive, loop[i] returns the
        index of the ith edge in the loop, and again the order matches the
        loop's orientation
    """
    def __init__( self, edges=None, orientation=0 ):
        """
        Creates a boundary loop from the given list of edges.

        If no edges are supplied, then creates an empty object with no data.
        In this case, one of the "set from" routines must be called on the
        boundary loop before performing any computations.

        If the optional orientation argument is not supplied, then the
        embedded loop will be assigned an arbitrary orientation. Otherwise,
        the supplied orientation must be one of the following:
        --> +1 if the first edge in the given list should be oriented from
            vertex 0 to vertex 1;
        --> -1 if the first edge should be oriented from vertex 1 to vertex 0;
            or
        --> 0 if this routine should be allowed to choose an arbitrary
            orientation.

        Raises NotLoop if the given list of edges does not form an embedded
        closed loop, or if the order of the edges in the given list does not
        match the order in which the edges appear in the loop.

        The given edges must all lie entirely on the boundary of the ambient
        triangulation; in other words, after constructing the boundary loop
        L, calling L.isBoundary() must return True. This condition is not
        checked, but some of the routines provided by this class might have
        undefined behaviour if this condition is not satisfied.

        Precondition:
        --> If supplied, the given list of edges must be nonempty, must
            consist of edges that all belong to the same 3-manifold
            triangulation T, and moreover all of these edges must lie
            entirely on the boundary of T.
        """
        super().__init__( edges, orientation )
        return

    def boundaryComponent(self):
        """
        Returns the boundary component containing this boundary loop.
        """
        return self._tri.edge( self[0] ).boundaryComponent()

    def shorten(self):
        """
        Shortens this boundary loop.

        In detail, if this boundary loop meets any *boundary* triangle F in
        exactly two distinct edges, then it can be shortened by replacing
        these two edges with the third edge of F. This routine performs such
        shortenings until no further shortening is possible.

        This routine raises BoundsDisc if self.triangulation() includes a
        boundary triangular face F that forms an embedded disc whose boundary
        is given by this boundary loop.

        Returns:
            True if and only if this boundary loop was successfully
            shortened. In the case where no shortening occurred, this
            boundary loop will remain entirely untouched.
        """
        # Can use the default implementation provided we supply an
        # implementation for _redirectCandidates().
        # Since _redirectCandidates() yields all boundary triangles incident
        # to this loop, this will raise BoundsDisc as promised in the
        # documentation.
        return super().shorten()

    def _redirectCandidates(self):
        """
        Yields candidate triangles of self.triangulation() across which the
        shorten() routine should attempt to redirect this loop.

        For a BoundaryLoop, every boundary triangle incident to the loop is a
        candidate. (We can only redirect along boundary triangles if we wish
        to ensure that the loop stays in the boundary.)
        """
        for ei in self:
            edge = self._tri.edge(ei)

            # Note that as a precondition, the current edge is assumed to be
            # a boundary edge. Thus, the incident boundary faces are at the
            # front and back of the current edge.
            front = edge.front()
            fFace = front.tetrahedron().triangle( front.vertices()[3] )
            yield fFace
            back = edge.back()
            bFace = back.tetrahedron().triangle( back.vertices()[2] )
            if bFace != fFace:
                yield bFace
        return
