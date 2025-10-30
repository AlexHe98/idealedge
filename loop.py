"""
Embedded loops in a 3-manifold triangulation, which play two main roles:
--> Ideal loops for representing torus boundary components of a 3-manifold.
--> Boundary loops on triangulation with real boundary.
"""
from regina import *
from moves import twoThree, threeTwo, twoZero, twoOne, fourFour
from insert import snapEdge, layerOn


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


#TODO Go through the entire class and its subclasses, and check what needs to
#   be modified to track orientations.
#TODO Find all uses of the "setFrom" routines, and see how the usage needs to
#   be modified to track orientations.
#TODO Find all uses of the IdealLoop() or BoundaryLoop() constructors, and
#   make sure that they track orientations (if necessary).
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
    def __init__( self, edges=None, orientation=None ):
        """
        Creates an embedded loop from the given list of edges.

        If no edges are supplied, then creates an empty object with no data.
        In this case, one of the "set from" routines must be called on the
        embedded loop before performing any computations.

        If the optional orientation argument is not supplied, then the
        embedded loop will be assigned an arbitrary orientation. Otherwise,
        the supplied orientation must be either +1 or -1; orientation +1
        means that the first edge e in the given list of edges (that is,
        e = edges[0]) is oriented from vertex 0 to vertex 1, whilst
        orientation -1 means that e is oriented from vertex 1 to vertex 0.

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

    def setFromEdges( self, edges, orientation=None ):
        """
        Sets this embedded loop using the given list of edges.

        If the optional orientation argument is not supplied, then the
        embedded loop will be assigned an arbitrary orientation. Otherwise,
        the supplied orientation must be either +1 or -1; orientation +1
        means that the first edge e in the given list of edges (that is,
        e = edges[0]) is oriented from vertex 0 to vertex 1, whilst
        orientation -1 means that e is oriented from vertex 1 to vertex 0.

        Raises NotLoop if the given list of edges does not form an embedded
        closed loop, or if the order of the edges in the given list does not
        match the order in which the edges appear in the loop.

        Precondition:
        --> The given list of edges is nonempty, and consists of edges that
            all belong to the same 3-manifold triangulation.
        """
        if ( ( orientation is not None ) and ( orientation != 1 ) and
            ( orientation != -1 ) ):
            raise ValueError(
                    "If supplied, orientation must be either +1 or -1." )
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

        The orientation of the cloned loop is guaranteed to match the
        orientation of the given loop.
        """
        if copyTri:
            tri = Triangulation3( loop.triangulation() )
        else:
            tri = loop.triangulation()
        edges = []
        for ei in loop:
            edges.append( tri.edge(ei) )
        self.setFromEdges( edges, loop.orientation() )
        return

    def setFromEdgeLocations( self, edgeLocations, orientation=None ):
        """
        Sets this embedded loop using the given list of edge locations.

        In detail, each edge location must be a pair (t, e), where t is a
        tetrahedron and e is an edge number from 0 to 5 (inclusive). Each
        tetrahedron must belong to the same triangulation.

        If the optional orientation argument is not supplied, then the
        embedded loop will be assigned an arbitrary orientation. Otherwise,
        the supplied orientation must be either +1 or -1; orientation +1
        means that the first edge e in the given list of edges (that is,
        e = edges[0]) is oriented from vertex 0 to vertex 1, whilst
        orientation -1 means that e is oriented from vertex 1 to vertex 0.

        Raises NotLoop if the given list of edges does not form an embedded
        closed loop, or if the order of the edges in the given list does not
        match the order in which the edges appear in the loop.

        Precondition:
        --> The given list of edge locations is nonempty.
        --> The tetrahedra given by the first entries of each edge location
            must all belong to the same 3-manifold triangulation.
        """
        edges = []
        for tet, edgeNum in edgeLocations:
            edges.append( tet.edge(edgeNum) )
        self.setFromEdges( edges, orientation )
        return

    #TODO Delete this entirely at a later date.
    def setFromLightweight( self, sig, edgeLocations ):
        """
        This routine is no longer available; use the new setFromBlueprint()
        routine instead.

        The old "lightweight description" consisted of an isomorphism
        signature and a list of edge locations. Unlike the new "picklable
        blueprint", this did not keep track of the orientation of the
        embedded loop, and only provided enough information to reconstruct
        the ambient triangulation up to combinatorial isomorphism.
        """
        raise NotImplementedError( "Use setFromBlueprint() instead." )

    def setFromBlueprint( self, sig, tetImages, facePerms,
                         edgeIndices, orientation ):
        """
        Sets this embedded loop using a picklable blueprint, as constructed
        by EmbeddedLoop.blueprint().
        """
        # Reconstruct the triangulation with the original labelling (not the
        # canonical labelling given by fromIsoSig()).
        tri = Triangulation3.fromIsoSig(sig)
        isom = Isomorphism3( tri.size() )
        for i in range( tri.size() ):
            isom.setSimpImage( i, tetImages[i] )
            isom.setFacePerm( i, Perm4( *facePerms[i] ) )

        # To directly clone the embedded loop, we need to recover the
        # original labelling.
        tri = isom.inverse()(tri)
        edges = [ tri.edge(ei) for ei in edgeIndices ]
        self.setFromEdges( edges, orientation )
        return

    #TODO What do we need to change to track orientations?
    def _setFromRenum( self, renum ):
        """
        Sets this embedded loop using the given edge renumbering map.

        This routine is for internal use only.
        """
        edges = []
        for ei in self:
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

    def vertexIndices(self):
        """
        Yields all indices of vertices incident to this embedded loop (in no
        particular order).
        """
        return iter( self._vertIndices )

    def clone(self):
        """
        Returns a clone of this embedded loop.

        The cloned loop will always be embedded in a copy of the
        triangulation containing this loop.
        """
        return EmbeddedLoop( self._cloneImpl() )

    def _cloneImpl(self):
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

    #TODO Delete this entirely at a later date.
    def lightweightDescription(self):
        """
        This routine is no longer available; use the new blueprint() routine
        instead.

        The old "lightweight description" consisted of an isomorphism
        signature and a list of edge locations. Unlike the new "picklable
        blueprint", this did not keep track of the orientation of the
        embedded loop, and only provided enough information to reconstruct
        the ambient triangulation up to combinatorial isomorphism.
        """
        raise NotImplementedError( "Use blueprint() instead." )

    def blueprint(self):
        """
        Returns a picklable blueprint for this embedded loop.

        In detail, this routine returns a tuple (S,T,F,E,O), where:
        --> S is the isomorphism signature for
            oldTri := self.triangulation();
        --> T is a list such that T[i] gives the index of the tetrahedron in
            newTri := Triangulation3.fromIsoSig(S) that corresponds to
            tetrahedron i of oldTri;
        --> F is a list such that F[i] is length-4 list with F[i][j] being
            the vertex number of newTri.tetrahedron( T[i] ) that corresponds
            to vertex j of oldTri.tetrahedron(i);
        --> E is the list of edge indices given by this embedded loop; and
        --> O is the orientation of this embedded loop.
        The returned blueprint can be used, via the setFromBlueprint()
        routine, to build a clone of this embedded loop.
        """
        sig, isom = self._tri.isoSigDetail()

        # Convert the isomorphism into the lists T and F.
        tetImages = []
        facePerms = []
        for i in range( self._tri.size() ):
            tetImages.append( isom.simpImage(i) )
            facePerms.append(
                    [ isom.facetPerm(i)[j] for j in range(4) ] )

        # For safety, return a copy of this loop's edge index list.
        return ( sig, tetImages, facePerms,
                list(self._edgeIndices), self.orientation() )

    def intersects( self, surf ):
        """
        Returns True if and only if this embedded loop has nonempty
        intersection with the given normal surface surf.

        Precondition:
        --> The given normal surface is embedded in self.triangulation().
        """
        for i in self:
            if surf.edgeWeight(i).safeLongValue() > 0:
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
            wt += surf.edgeWeight(i).safeLongValue()
        return wt

    #TODO Delete this entirely at a later date.
    def components( self, surf ):
        """
        This routine is no longer available; use the new splitArcs() routine
        instead.

        Apart from the improvement in the name, the new splitArcs() routine
        also keeps track of the orientation of the embedded loop.
        """
        raise NotImplementedError()

    def splitArcs( self, surf ):
        """
        Returns a list describing the arcs into which the given normal
        surface surf splits this embedded loop.

        In detail, each item of the returned list is a list of edge segments
        that all belong to the same arc, satisfying the following conditions:
        --> The edge segments appear in the same order as they do when we
            traverse this embedded loop.
        --> Each edge segment in L is encoded as a triple (ei, n, o) such
            that:
            --- ei is the index of the edge containing the segment in
                question;
            --- n is the segment number (from 0 to w, inclusive, where w is
                the weight of the given surface on edge ei); and
            --- o is +1 if edge ei is oriented from vertex 0 to vertex 1, and
                -1 if edge ei is oriented from vertex 1 to vertex 0.
        Note that for each edge e, the segments are numbered in ascending
        order from the segment incident to e.vertex(0) to the segment
        incident to e.vertex(1).

        Note also that the order of the arcs need not be the same as the
        order in which they appear as we traverse this loop. Only the order
        of edge segments within each individual arc is guaranteed to match
        the traversal order.
        """
        # We find all the arcs by simply walking around the loop. Take the
        # first arc to be the one that begins *after* the first point at
        # which this loop gets split by the given surf. Thus, our walk starts
        # in the middle of the last arc, so we need to make sure to remember
        # all the segments of the last arc.
        lastArc = []
        splitIndex = None
        for i in range( len(self) ):
            edgeIndex = self._edgeIndices[i]
            orientation = self.edgeOrientation(i)
            wt = surf.edgeWeight(edgeIndex).safeLongValue()
            if wt > 0:
                # We found the point at which the first arc begins.
                if self._tails[i] == 0:
                    lastArc.append( ( edgeIndex, 0, orientation ) )
                    headSeg = ( edgeIndex, wt, orientation )
                else:
                    lastArc.append( ( edgeIndex, wt, orientation ) )
                    headSeg = ( edgeIndex, 0, orientation )
                splitIndex = i
                break
            else:
                # We are still in the middle of the last arc.
                lastArc.append( ( edgeIndex, 0, orientation ) )
        if splitIndex is None:
            # If this loop is disjoint from the surface, then there is only
            # one arc, and we have already found all the constituent segments
            # of this arc.
            return [lastArc]

        # The given surf splits this embedded loop into multiple arcs, so we
        # need to do a bit more work.
        arcs = []
        while splitIndex is not None:
            orientation = headSeg[-1]
            for seg in range( 1, wt ):
                # For wt >= 2, we get a sequence of short arcs given by
                # type-2 segments.
                #
                # If orientation == -1 and wt > 2, then this will add new
                # arcs in the "wrong" order.
                arcs.append( [ ( edgeIndex, seg, orientation ) ] )

            # We now need to find all the segments that comprise the next
            # (long) arc.
            nextArc = [headSeg]
            continuation = splitIndex + 1

            # Unless we have already returned to the last arc, we must
            # eventually find another split point at which the next arc
            # begins.
            splitIndex = None
            for i in range( continuation, len(self) ):
                edgeIndex = self._edgeIndices[i]
                orientation = self.edgeOrientation(i)
                wt = surf.edgeWeight(edgeIndex).safeLongValue()
                if wt > 0:
                    # Found the next split point.
                    if self._tails[i] == 0:
                        nextArc.append( ( edgeIndex, 0, orientation ) )
                        headSeg = ( edgeIndex, wt, orientation )
                    else:
                        nextArc.append( ( edgeIndex, wt, orientation ) )
                        headSeg = ( edgeIndex, 0, orientation )
                    splitIndex = i
                    arcs.append(nextArc)
                    break
                else:
                    # We are still in the middle of the current arc.
                    nextArc.append( ( edgeIndex, 0, orientation ) )

        # Don't forget to include the last arc.
        arcs.append( [ *nextArc, *lastArc ] )
        return arcs

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

    #TODO Check what needs to be done for orientations for everything below
    #   this point.
    #TODO WORKING HERE.

    def _shortenImpl(self):
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

        This routine might raise BoundsDisc.

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
        _shortenImpl() routine should attempt to redirect this loop.

        The EmbeddedLoop base class does not implement this routine, so
        subclasses that require this routine must provide an implementation.
        """
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
        newEdges = [ self._tri.edge(i) for i in self ]
        self.setFromEdges(newEdges)
        return True

    def _minimiseBoundaryImpl(self):
        """
        Ensures that the triangulation containing this embedded loop has the
        smallest possible number of boundary triangles, potentially adding
        tetrahedra to do this.

        The default implementation of this routine requires the following
        helper routines, which are *not* fully implemented by default:
        --> _shortenImpl()
        --> _findBoundaryMove()
        Thus, subclasses that require this routine must either:
        --> override this routine; or
        --> supply implementations for the aforementioned helper routines.
        In the latter case, see the documentation for each respective helper
        routine for details on the behaviour that must be implemented.

        A side-effect of calling this routine is that it will shorten this
        embedded loop if possible.

        This routine might raise BoundsDisc.

        The following are guaranteed to hold once this routine is finished:
        --> Every 2-sphere boundary component will have exactly two triangles
            and three vertices.
        --> Every projective plane boundary component will have exactly two
            triangles and two vertices.
        --> Every other boundary component will have exactly one vertex.

        The changes that this routine performs can always be expressed using
        only the following operations:
        --> Shortening this loop by redirecting it across triangular faces.
        --> Close book moves and/or layerings on self.triangulation().
        In particular, this routine never creates new vertices, and it never
        creates a non-vertex-linking normal disc or 2-sphere if there was not
        one before.

        If the triangulation containing this loop is currently oriented, then
        this routine guarantees to preserve the orientation.

        Adapted from Regina's Triangulation3.minimiseBoundary().

        Precondition:
        --> The ambient triangulation (i.e., self.triangulation()) is valid.

        Returns:
            True if and only if this loop or its ambient triangulation were
            changed. In other words, a return value of False indicates that:
            (1) this loop could not be shortened; and
            (2) every boundary component of the ambient triangulation was
                already minimal to begin with.
        """
        # Simplify the boundary of self._tri to use as few triangles as
        # possible.
        changed = False
        while True:
            # Shorten this loop to minimise the number of special cases.
            if self._shortenImpl():     # Might raise BoundsDisc.
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
                edge = layerOn(edge).edge(5)
            self._tri.closeBook( edge, False, True )
            self.setFromEdgeLocations(edgeLocations)
        return

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

        If the triangulation containing this loop is currently oriented, then
        this routine guarantees to preserve the orientation.

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

        If the triangulation containing this loop is currently oriented, then
        this routine guarantees to preserve the orientation.

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

        If the triangulation containing this loop is currently oriented, then
        this routine guarantees to preserve the orientation.

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


#TODO Update class documentation to mention tracking of orientation.
class IdealLoop(EmbeddedLoop):
    """
    A sequence of edges representing an embedded ideal loop in the interior
    of a 3-manifold triangulation.

    Some of the routines provided by this class might fail if the ideal loop
    bounds an embedded disc in the ambient triangulation (though these
    routines might nevertheless succeed in spite of the existence of such a
    disc). This class raises BoundsDisc whenever such a failure occurs.

    A core feature of this class is that it effectively stores a list of edge
    *indices* corresponding to the edges of the ideal loop. Thus, for any
    instance loop of this class, the following functionality is available:
    --> (e in loop) is True if and only if loop.triangulation().edge(e) is an
        edge in the loop
    --> len(loop) is the number of edges in the loop
    --> for i between 0 and (len(loop) - 1), inclusive, loop[i] returns the
        index of the ith edge in the loop
    --> iterating through the loop yields all the edge indices in order
    """
    def __init__( self, edges=None, orientation=None ):
        """
        Creates an ideal loop from the given list of edges.

        If no edges are supplied, then creates an empty object with no data.
        In this case, one of the "set from" routines must be called on the
        ideal loop before performing any computations.

        If the optional orientation argument is not supplied, then the
        embedded loop will be assigned an arbitrary orientation. Otherwise,
        the supplied orientation must be either +1 or -1; orientation +1
        means that the first edge e in the given list of edges (that is,
        e = edges[0]) is oriented from vertex 0 to vertex 1, whilst
        orientation -1 means that e is oriented from vertex 1 to vertex 0.

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
        for ei in self:
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

        This routine might raise BoundsDisc.

        Returns:
            True if and only if this ideal loop was successfully shortened.
            In the case where no shortening occurred, this ideal loop will
            remain entirely untouched.
            True if and only if this ideal loop was successfully shortened.
            Otherwise, this ideal loop will not be modified at all.
        """
        # Can use the default implementation provided we supply an
        # implementation for _redirectCandidates().
        return self._shortenImpl()  # Might raise BoundsDisc.

    def _redirectCandidates(self):
        # Any triangle that is incident to this ideal loop is a candidate to
        # use for redirecting.
        for ei in self:
            edge = self._tri.edge(ei)

            # Yield *all* triangles incident to current edge.
            for emb in edge.embeddings():
                yield emb.tetrahedron().triangle( emb.vertices()[3] )
            if edge.isBoundary():
                yield emb.tetrahedron().triangle( emb.vertices()[2] )
        return

    def minimiseBoundary(self):
        """
        Ensures that the triangulation containing this ideal loop has the
        smallest possible number of boundary triangles, potentially adding
        tetrahedra to do this.

        A side-effect of calling this routine is that it will shorten this
        ideal loop if possible.

        This routine might raise BoundsDisc.

        The following are guaranteed to hold once this routine is finished:
        --> Every 2-sphere boundary component will have exactly two triangles
            and three vertices.
        --> Every projective plane boundary component will have exactly two
            triangles and two vertices.
        --> Every other boundary component will have exactly one vertex.

        The changes that this routine performs can always be expressed using
        only the following operations:
        --> Shortening this loop by redirecting it across triangular faces.
        --> Close book moves and/or layerings on self.triangulation().
        In particular, this routine never creates new vertices, and it never
        creates a non-vertex-linking normal disc or 2-sphere if there was not
        one before.

        If the triangulation containing this loop is currently oriented, then
        this routine guarantees to preserve the orientation.

        Adapted from Regina's Triangulation3.minimiseBoundary().

        Precondition:
        --> The ambient triangulation (i.e., self.triangulation()) is valid.

        Returns:
            True if and only if this loop or its ambient triangulation were
            changed. In other words, a return value of False indicates that:
            (1) this loop could not be shortened; and
            (2) every boundary component of the ambient triangulation was
                already minimal to begin with.
        """
        # Can use the default implementation provided we supply an
        # implementation for _findBoundaryMove().
        return self._minimiseBoundaryImpl()

    def _findBoundaryMove(self):
        # Precondition:
        #   --> This loop cannot be shortened.

        # Find a boundary component that is not yet minimal.
        for bc in self._tri.boundaryComponents():
            if bc.countTriangles() <= 2 or bc.countVertices() <= 1:
                # Current boundary component is already minimal.
                continue

            # First try to find a close book move, which does not increase
            # the number of tetrahedra.
            for e in bc.edges():
                if self._tri.closeBook( edge, True, False ):
                    return ( edge,
                            False,  # Close book without layering.
                            self._edgeIndices )

            # We could not find a close book move.
            # In this case, because bc is non-minimal, there must be a
            # boundary edge e that joins two distinct vertices, and we can
            # simplify bc by layering across e and then performing a close
            # book move on the newly layered edge.
            for e in bc.edges():
                if edge.vertex(0) == edge.vertex(1):
                    continue

                # The layering is illegal if this edge is incident to the
                # same boundary triangle F on both sides (rather than two
                # distinct triangles). But in that scenario, F forms a disc,
                # and there must be a close book move available on the edge b
                # that forms the boundary of this disc. Thus, when we reach
                # this point in the code, we can guarantee that the layering
                # is legal.
                return ( edge,
                        True,   # Layer before performing close book.
                        self._edgeIndices )

            # We should never reach this point.
            raise RuntimeError(
                    "_findBoundaryMove() failed unexpectedly." )

        # If we fell out of the boundary component loop, then all boundary
        # components are minimal.
        return None

    def minimiseVertices(self):
        """
        Ensures that the triangulation containing this ideal loop has the
        smallest possible number of vertices for the 3-manifold that it
        represents, potentially adding tetrahedra to do this.

        A side-effect of calling this routine is that it will shorten this
        ideal loop if possible.

        This routine might raise BoundsDisc.

        The following are guaranteed to hold once this routine is finished:
        --> If the ambient triangulation is closed, then it will have
            precisely one vertex.
        --> If the ambient triangulation has real boundary, then:
            --- there will be exactly one internal vertex;
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

        If the triangulation containing this loop is currently oriented, then
        this routine guarantees to preserve the orientation.

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
        # Can use the default implementation provided we supply an
        # implementation for _findSnapEdge().
        return self._minimiseVerticesImpl()

    def _findSnapEdge(self):
        # Precondition:
        #   --> This loop cannot be shortened.
        #   --> The boundary of self._tri has been minimised.

        # Find a suitable edge on which to perform a snap edge move. We
        # minimise the number of special cases by prioritising edges
        # belonging to the ideal loop.
        if len(self) > 1:
            # This loop is supposed to lie entirely in the interior of the
            # ambient triangulation, so it should be legal to snap any edge
            # of this loop; here, we choose the last edge in the loop.
            #
            # It is safe to directly modify self._edgeIndices since this will
            # need to be updated anyway.
            lastEdgeIndex = self._edgeIndices.pop()
            return ( self._tri.edge(lastEdgeIndex), self._edgeIndices )
        else:
            # At this point, this loop is one-vertex and the boundary has
            # been minimised, but there might still be some other internal
            # vertices that we can remove.
            for edge in self._tri.edges():
                if not snapEdge( edge, True, False ):
                    # Snap edge is not legal.
                    continue

                # The snap edge move is legal, but we only want to perform
                # this move if it will remove an internal vertex that is not
                # incident to this ideal loop.
                for i in range(2):
                    v = edge.vertex(i)
                    if v.isBoundary() or v.index() in self._vertIndices:
                        continue

                    # We can use this snap edge move to remove v!
                    return ( edge, self._edgeIndices )
        return

    def simplifyBasic(self):
        """
        Uses 2-0 edge and 2-1 edge moves to monotonically reduce the number
        of tetrahedra in the ambient triangulation, while leaving this ideal
        loop untouched.

        There should usually be no need to call this routine directly, since
        the functionality is subsumed by the more powerful simplify() and
        simplifyMonotonic() routines.

        This routine might raise BoundsDisc.

        If the triangulation containing this loop is currently oriented, then
        this routine guarantees to preserve the orientation.

        Adapted from SnapPea's check_for_cancellation().

        Returns:
            True if and only if the ambient triangulation was successfully
            simplified. Otherwise, the ambient triangulation will not be
            modified at all.
        """
        # Do not include 3-2 moves.
        # Might raise BoundsDisc.
        return self._simplifyMonotonicImpl(False)

    def simplifyMonotonic(self):
        """
        Uses 2-0 edge, 2-1 edge, and 3-2 moves to monotonically reduce the
        number of tetrahedra in the ambient triangulation, while leaving this
        ideal loop untouched.

        There should usually be no need to call this routine directly, since
        the functionality is subsumed by the more powerful simplify()
        routine.

        This routine might raise BoundsDisc.

        If the triangulation containing this loop is currently oriented, then
        this routine guarantees to preserve the orientation.

        Adapted from Regina's Triangulation3.simplifyToLocalMinimum().

        Returns:
            True if and only if the ambient triangulation was successfully
            simplified. Otherwise, the ambient triangulation will not be
            modified at all.
        """
        # Include 3-2 moves.
        # Might raise BoundsDisc.
        return self._simplifyMonotonicImpl(True)

    def simplify(self):
        """
        Attempts to simplify this ideal loop.

        This routine uses minimiseVertices() and simplifyMonotonic(), in
        combination with random 4-4 moves that leave this loop untouched.
        
        Although this routine works very well most of the time, it can
        occasionally get stuck in a "well" that can only be escaped by
        increasing the number of tetrahedra. In such cases, it might be
        useful to try to escape using the randomise() routine.

        This routine might raise BoundsDisc.

        If the triangulation containing this loop is currently oriented, then
        this routine guarantees to preserve the orientation.

        Adapted from Regina's Triangulation3.intelligentSimplify().

        Warning:
            --> Running this routine multiple times on the same loop may
                return different results, since the implementation makes
                random decisions.

        Returns:
            True if and only if this loop was successfully simplified.
            Otherwise, this loop will not be modified at all.
        """
        return self._simplifyImpl()

    def randomise(self):
        """
        Attempts to randomly retriangulate this ideal loop.

        This routine works by performing lots of random 2-3 moves, before
        attempting to simplify the ambient triangulation again. As such, it
        is often useful for escaping "wells" when the simplify() routine gets
        stuck.

        This routine might raise BoundsDisc.

        If the triangulation containing this loop is currently oriented, then
        this routine guarantees to preserve the orientation.

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


#TODO Update class documentation to mention tracking of orientation.
class BoundaryLoop(EmbeddedLoop):
    """
    A sequence of edges representing an embedded loop on the boundary of a
    3-manifold triangulation.

    Some of the routines provided by this class might fail if the boundary
    loop bounds an embedded disc in the ambient triangulation (though these
    routines might nevertheless succeed in spite of the existence of such a
    disc). This class raises BoundsDisc whenever such a failure occurs.

    A core feature of this class is that it effectively stores a list of edge
    *indices* corresponding to the edges of the boundary loop. Thus, for any
    instance loop of this class, the following functionality is available:
    --> (e in loop) is True if and only if loop.triangulation().edge(e) is an
        edge in the loop
    --> len(loop) is the number of edges in the loop
    --> for i between 0 and (len(loop) - 1), inclusive, loop[i] returns the
        index of the ith edge in the loop
    --> iterating through the loop yields all the edge indices in order
    """
    def __init__( self, edges=None, orientation=None ):
        """
        Creates a boundary loop from the given list of edges.

        If no edges are supplied, then creates an empty object with no data.
        In this case, one of the "set from" routines must be called on the
        boundary loop before performing any computations.

        If the optional orientation argument is not supplied, then the
        embedded loop will be assigned an arbitrary orientation. Otherwise,
        the supplied orientation must be either +1 or -1; orientation +1
        means that the first edge e in the given list of edges (that is,
        e = edges[0]) is oriented from vertex 0 to vertex 1, whilst
        orientation -1 means that e is oriented from vertex 1 to vertex 0.

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

    def clone(self):
        """
        Returns a clone of this boundary loop.

        The cloned boundary loop will always be embedded in a copy of the
        triangulation containing this boundary loop.
        """
        return BoundaryLoop( self._cloneImpl() )

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

        This routine might raise BoundsDisc.

        Returns:
            True if and only if this boundary loop was successfully
            shortened. In the case where no shortening occurred, this
            boundary loop will remain entirely untouched.
        """
        # Can use the default implementation provided we supply an
        # implementation for _redirectCandidates().
        return self._shortenImpl()  # Might raise BoundsDisc.

    def _redirectCandidates(self):
        # For boundary loops, we only want to redirect along *boundary*
        # triangles (to ensure that the loop stays in the boundary).
        for ei in self:
            edge = self._tri.edge(ei)

            # Note that as a precondition, the current edge is assumed to be
            # a boundary edge. Thus, the incident boundary faces are at the
            # front and back of the current edge.
            front = edge.front()
            yield front.tetrahedron().triangle( front.vertices()[3] )
            back = edge.back()
            yield back.tetrahedron().triangle( back.vertices()[2] )
        return

    def minimiseBoundary(self):
        """
        Ensures that the triangulation containing this boundary loop has the
        smallest possible number of boundary triangles, potentially adding
        tetrahedra to do this.

        A side-effect of calling this routine is that it will shorten this
        boundary loop if possible.

        This routine might raise BoundsDisc.

        The following are guaranteed to hold once this routine is finished:
        --> Every 2-sphere boundary component will have exactly two triangles
            and three vertices.
        --> Every projective plane boundary component will have exactly two
            triangles and two vertices.
        --> Every other boundary component will have exactly one vertex.

        The changes that this routine performs can always be expressed using
        only the following operations:
        --> Shortening this loop by redirecting it across triangular faces.
        --> Close book moves and/or layerings on self.triangulation().
        In particular, this routine never creates new vertices, and it never
        creates a non-vertex-linking normal disc or 2-sphere if there was not
        one before.

        If the triangulation containing this loop is currently oriented, then
        this routine guarantees to preserve the orientation.

        Adapted from Regina's Triangulation3.minimiseBoundary().

        Precondition:
        --> The ambient triangulation (i.e., self.triangulation()) is valid.

        Returns:
            True if and only if this loop or its ambient triangulation were
            changed. In other words, a return value of False indicates that:
            (1) this loop could not be shortened; and
            (2) every boundary component of the ambient triangulation was
                already minimal to begin with.
        """
        # Can use the default implementation provided we supply an
        # implementation for _findBoundaryMove().
        return self._minimiseBoundaryImpl()

    def _findBoundaryMove(self):
        # Exceptions:
        #   --> Might raise BoundsDisc.
        #
        # Precondition:
        #   --> This loop cannot be shortened.

        # Prioritise moves that reduce the length of this ideal loop. If
        # possible, use close book moves so that we do not introduce too many
        # new tetrahedra.
        if len(self) > 1 and self.boundaryComponent().countTriangles() > 2:
            # Try to find a close book move.
            if len(self) > 2:
                for edge in self.boundaryComponent().edges():
                    # Check eligibility of close book move, but do *not*
                    # perform yet.
                    if not self._tri.closeBook( edge, True, False ):
                        continue

                    # Does this close book move reduce the length?
                    ftet = edge.front().tetrahedron()
                    fver = edge.front().vertices()
                    btet = edge.back().tetrahedron()
                    bver = edge.back().vertices()
                    for v in range(2):
                        fei = ftet.edge( fver[v], fver[2] ).index()
                        bei = btet.edge( bver[v], bver[3] ).index()
                        if ( fei in self ) and ( bei in self ):
                            # Since this loop cannot be shortened, the
                            # indices fei and bei correspond to the *only*
                            # edges of this loop that are incident to the
                            # triangles involved in the close book move.
                            # After performing the move, we can remove these
                            # two edges from the loop.
                            newEdgeIndices = [ ei for ei in self
                                    if ei not in { fei, bei } ]
                            return ( edge,
                                    False,  # Close book without layering.
                                    newEdgeIndices )

            # Resort to layering a new tetrahedron to facilitate a close book
            # move that effectively removes one of the edges from this loop.
            # This operation is guaranteed to be legal for any edge belonging
            # to this loop; here we choose to perform it on the last edge.
            #
            # It is safe to directly modify self._edgeIndices since this will
            # need to be updated anyway.
            lastEdgeIndex = self._edgeIndices.pop()
            return ( self._tri.edge(lastEdgeIndex),
                    True,   # Layer before performing close book.
                    self._edgeIndices )

        # At this point, if the boundary component containing the loop is not
        # yet minimal, then we at least know that the loop consists only of a
        # single edge e. Our goal now is to minimise boundary components
        # without touching e. 
        for bc in self._tri.boundaryComponents():
            if bc.countTriangles() <= 2 or bc.countVertices() <= 1:
                continue

            # First try to find a close book move.
            for edge in bc.edges():
                if edge.index() == self[0]:
                    # Leave the loop untouched.
                    continue
                if self._tri.closeBook( edge, True, False ):
                    return ( edge,
                            False,  # Close book without layering.
                            self._edgeIndices )

            # We could not find a suitable close book move, so our plan now
            # is to find a boundary edge e that joins two distinct vertices.
            # We can layer over such an edge e, and then simplify the
            # boundary using a close book move on the newly layered edge.
            for edge in bc.edges():
                if edge.vertex(0) == edge.vertex(1):
                    continue

                # The layering is illegal if this edge is incident to the
                # same boundary triangle F on both sides (rather than two
                # distinct triangles). But in that scenario, F forms a disc,
                # and there must be a close book move available on the edge b
                # that forms the boundary of this disc. The only possible
                # close book that we haven't already ruled out is the one we
                # were trying to avoid; in other words, b must coincide with
                # the loop that we are trying to preserve.
                front = edge.front()
                back = edge.back()
                if ( front.tetrahedron().triangle( front.vertices()[3] ) ==
                        back.tetrahedron().triangle( back.vertices()[2] ) ):
                    raise BoundsDisc()
                else:
                    return ( edge,
                            True,   # Layer before performing close book.
                            self._edgeIndices )

            # We should never reach this point.
            raise RuntimeError(
                    "_findBoundaryMove() failed unexpectedly." )

        # If we fell out of the boundary component loop, then all boundary
        # components are minimal.
        return None

    def minimiseVertices(self):
        """
        Ensures that the triangulation containing this boundary loop has the
        smallest possible number of vertices for the 3-manifold that it
        represents, potentially adding tetrahedra to do this.

        A side-effect of calling this routine is that it will shorten this
        boundary loop if possible.

        This routine might raise BoundsDisc.

        The following are guaranteed to hold once this routine is finished:
        --> If the ambient triangulation is closed, then it will have
            precisely one vertex.
        --> If the ambient triangulation has real boundary, then:
            --- there will be no internal vertices;
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

        If the triangulation containing this loop is currently oriented, then
        this routine guarantees to preserve the orientation.

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
        # Can use the default implementation provided we supply an
        # implementation for _findSnapEdge().
        return self._minimiseVerticesImpl()

    def _findSnapEdge(self):
        # Precondition:
        #   --> This loop cannot be shortened.
        #   --> The boundary of self._tri has been minimised.

        # Find a suitable edge on which to perform a snap edge move (just
        # check whether the move is legal, do not perform yet).
        for edge in self._tri.edges():
            if snapEdge( edge, True, False ):
                # This loop should lie entirely in the boundary, so the snap
                # edge move should not change the loop.
                return ( edge, self._edgeIndices )
        return

    def simplifyMonotonic(self):
        """
        Uses 2-0 edge, 2-1 edge, and 3-2 moves to monotonically reduce the
        number of tetrahedra in the ambient triangulation, while leaving this
        boundary loop untouched.

        There should usually be no need to call this routine directly, since
        the functionality is subsumed by the more powerful simplify()
        routine.

        This routine might raise BoundsDisc.

        If the triangulation containing this loop is currently oriented, then
        this routine guarantees to preserve the orientation.

        Adapted from Regina's Triangulation3.simplifyToLocalMinimum().

        Returns:
            True if and only if the ambient triangulation was successfully
            simplified. Otherwise, the ambient triangulation will not be
            modified at all.
        """
        # Include 3-2 moves.
        # Might raise BoundsDisc.
        return self._simplifyMonotonicImpl(True)

    def simplify(self):
        """
        Attempts to simplify this boundary loop.

        This routine uses minimiseVertices() and simplifyMonotonic(), in
        combination with random 4-4 moves (which leave this loop untouched).

        This routine might raise BoundsDisc.

        If the triangulation containing this loop is currently oriented, then
        this routine guarantees to preserve the orientation.

        Adapted from Regina's Triangulation3.intelligentSimplify().

        Warning:
            --> Running this routine multiple times on the same loop may
                return different results, since the implementation makes
                random decisions.

        Returns:
            True if and only if this loop was successfully simplified.
            Otherwise, this loop will not be modified at all.
        """
        return self._simplifyImpl()
