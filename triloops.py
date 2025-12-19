"""
3-manifold triangulations containing disjoint unions embedded loops.

The classes here provide methods for simplifying the ambient 3-manifold
triangulation, while preserving the topological embedding of the loops.
"""
from regina import *
from moves import twoThree, threeTwo, twoZero, twoOne, fourFour
from insert import snapEdge, layerOn
from loop import EmbeddedLoop, IdealLoop, BoundaryLoop
from loopaux import embeddingsFromEdgeIndices


class TriangulationWithEmbeddedLoops:
    """
    A 3-manifold triangulation containing a disjoint union of EmbeddedLoop
    objects.

    This is a base class that implements common functionality for the
    EdgeIdealTriangulation and TriangulationWithBoundaryLoops classes.
    Although this base class can be instantiated, the functionality it offers
    is much less complete than its aforementioned subclasses.

    This class has two core features:
    (1) It provides methods to simplify the ambient 3-manifold triangulation,
        while ensuring that the topological embedding of the union of loops is
        always preserved.
    (2) It acts as a container of EmbeddedLoop objects, which are indexed in
        an arbitrary order (but the order is kept consistent no matter how
        much the ambient triangulation is simplified). In detail, for any
        instance triLoops of this class:
        --> (el in triLoops) is True if and only if el is one of the
            EmbeddedLoop objects that is embedded in the ambient triangulation
            (note that equality of EmbeddedLoop objects is determined by their
            location in memory, and so for instance clones will not be
            considered equal)
        --> len(triLoops) is the number of EmbeddedLoop objects that are
            embedded in the ambient triangulation
        --> iterating through triLoops yields all the EmbeddedLoop objects
            that are embedded in the ambient triangulation, in the order in
            which they are indexed
        --> for i between 0 and (len(triLoops) - 1), inclusive, triLoops[i]
            returns the EmbeddedLoop at index i
    """
    # Base class for loops embedded in the ambient triangulation.
    _LOOP_CLASS = EmbeddedLoop

    def __init__( self, loops ):
        """
        Creates a disjoint union of the given collection of loops.

        Precondition:
        --> loops is nonempty.
        --> The elements of loops are all EmbeddedLoop objects lying
            disjointly inside the same ambient 3-manifold triangulation.
        """
        self.setFromLoops(loops)
        return

    def setFromLoops( self, loops ):
        """
        Sets this to be a triangulation containing a disjoint union of the
        given loops.

        Precondition:
        --> loops is nonempty.
        --> The elements of loops are all EmbeddedLoop objects lying
            disjointly inside the same ambient 3-manifold triangulation.
        """
        self._loops = list(loops)
        self._tri = self._loops[0].triangulation()
        return

    def clone(self):
        """
        Returns a clone of this triangulation with embedded loops.

        The clone will always use a copy of self.triangulation()
        """
        # We use the built-in type() function to make sure that subclasses
        # will construct clones of the correct type.
        cloneTri = Triangulation3(self._tri)
        cloneLoops = []
        for loop in self._loops:
            cloneLoops.append( type(loop)(
                loop._cloneImpl(cloneTri), loop.orientation() ) )
        return type(self)(cloneLoops)

    @classmethod
    def fromBlueprint( cls, triEncoding, *loops ):
        """
        Constructs a triangulation with embedded loops using a picklable
        blueprint, as constructed by either EmbeddedLoop.blueprint() or
        TriangulationWithEmbeddedLoops.blueprint().
        """
        tri = Triangulation3.tightDecoding(triEncoding)
        embLoops = []
        i = 0
        if len(loops) % 2 == 1:
            raise ValueError("Invalid blueprint")
        while i < len(loops):
            edgeIndices = loops[i]
            orientation = loops[i+1]
            i += 2
            edges = [ tri.edge(ei) for ei in edgeIndices ]
            embLoops.append( cls._LOOP_CLASS( edges, orientation ) )
        return cls(embLoops)

    def setFromEdgeEmbeddings( self, data ):
        """
        Sets this triangulation with embedded loops using the given data.

        In detail, each item in data is a pair (E, D) that specifies an
        EmbeddedLoop in a triangulation T as follows:
        --> E is a list of edge embeddings.
        --> D specifies the orientation of the EmbeddedLoop: 
            --- It is +1 if the first edge of the loop should be oriented from
                vertex 0 to vertex 1 (here, vertex numbers are with respect to
                the edge embedding, which might differ from the vertex numbers
                of the underlying edge if the ambient triangulation has been
                modified since the edge embedding was constructed).
            --- It is -1 if the first edge of the loop should be oriented from
                vertex 1 to vertex 0.
            --- It is 0 if this routine is allowed to choose an arbitrary
                orientation on the loop.
        All of the edge embeddings must reference tetrahedra in the same
        3-manifold triangulation.

        Raises NotLoop if any of the given lists of edge embeddings does not
        describe an embedded closed loop, or if the order of the edges in the
        list does not match the order in which the edges appear in the loop.

        Precondition:
        --> Each given list of edge embeddings is nonempty.
        --> All the edge embeddings must reference tetrahedra belonging to the
            same 3-manifold triangulation.
        """
        embLoops = []
        for edgeEmbeddings, orientation in data:
            embLoops.append(
                    self._LOOP_CLASS.fromEdgeEmbeddings(
                        edgeEmbeddings, orientation ) )
        self.setFromLoops(embLoops)
        return

    def _edgeEmbeddingsData( self, *, remove=set() ):
        """
        Returns data to reconstruct the embedded loops using the
        setFromEdgeEmbeddings() routine.

        See setFromEdgeEmbeddings() for more details on the data structure.

        The optional remove parameter allows for some edge indices to be
        removed from some of the loops. This is useful if we wish to modify
        the ambient triangulation in a way that would remove the corresponding
        edges from the loops; in such a situation, the returned data can be
        used to correctly reconstruct the loops in the new triangulation that
        results from performing the move.
        """
        ans = []
        for embLoop in self:
            # Get the data for the current embLoop, with edges removed if
            # necessary.
            edgeIndices = []
            orientation = None
            for ind, edgeInd in enumerate(embLoop):
                if edgeInd in remove:
                    continue
                edgeIndices.append(edgeInd)
                if orientation is None:
                    # Found the first non-removed edge. This is the one that
                    # will determine the new orientation.
                    orientation = embLoop.edgeOrientation(ind)

            # Add the data for the current embLoop to ans.
            embeddings = embeddingsFromEdgeIndices(
                    self._tri, edgeIndices )
            ans.append( ( embeddings, orientation ) )
        return ans

    def __len__(self):
        return len(self._loops)

    def __contains__( self, embLoop ):
        return embLoop in self._loops

    def __iter__(self):
        return iter( self._loops )

    def __getitem__( self, index ):
        return self._loops[index]

    def triangulation(self):
        """
        Returns the ambient triangulation.
        """
        return self._tri

    def isBoundary(self):
        """
        Does the union of embedded loops lie entirely in the boundary of the
        ambient triangulation?
        """
        for embLoop in self:
            if not embLoop.isBoundary():
                return False
        return True

    def isInternal(self):
        """
        Does the union of embedded loops lie entirely in the interior of the
        ambient triangulation?
        """
        for embLoop in self:
            if not embLoop.isInternal():
                return False
        return True

    def countLoopEdges(self):
        """
        Counts the number of edges in the union of embedded loops.
        """
        total = 0
        for embLoop in self:
            total += len(embLoop)
        return total

    def loopEdgeIndices(self):
        """
        Returns the set of all edge indices involved in the embedded loops.
        """
        ans = set()
        for embLoop in self:
            ans = ans.union(embLoop)
        return ans

    def loopVertexIndices(self):
        """
        Returns the set of all vertex indices involved in the embedded loops.
        """
        ans = set()
        for embLoop in self:
            ans = ans.union( embLoop.vertexIndices() )
        return ans

    def blueprint(self):
        """
        Returns a picklable blueprint for this triangulation with embedded
        loops.

        In detail, this routine returns a tuple
            (T, E<0>, D<0>, ..., E<L-1>, D<L-1>),
        where:
        --> L = len(self).
        --> T is Regina's tight encoding of self.triangulation().
        --> E<i> is (a copy of) the list of edge indices given by the embedded
            loop at index i.
        --> D<i> is the orientation of the embedded loop at index i.
        The returned blueprint can be used, via the class method
        fromBlueprint(), to build a clone of this triangulation with embedded
        loops.
        """
        ans = [ self._tri.tightEncoding() ]
        for loop in self:
            ans.append( loop.edgeIndices() )
            ans.append( loop.orientation() )
        return tuple(ans)

    #TODO Do we need intersects() and weight()?

    def intersects( self, surf ):
        """
        Returns True if and only if the union of embedded loops has nonempty
        intersection with the given normal surface surf.

        Precondition:
        --> The given normal surface is embedded in self.triangulation().
        """
        for embLoop in self:
            if embLoop.intersects(surf):
                return True
        return False

    def weight( self, surf ):
        """
        Returns the number of times the union of embedded loops intersects
        the given normal surface surf.

        Precondition:
        --> The given normal surface is embedded in self.triangulation().
        """
        wt = 0
        for embLoop in self:
            wt += embLoop.weight()
        return wt

    #TODO Need a version of splitArcs() that returns enough information to
    #   track all of the effects of crushing that we care about.

    #TODO Implement simplification, and remove the old simplification code
    #   from EmbeddedLoop.

    def shorten(self):
        """
        Shortens the union of embedded loops by looking for triangles that
        intersect a loop in two edges, and redirecting this loop to use the
        third edge.

        If no shortening is possible, then this triangulation with embedded
        loops will remain entirely untouched.

        If one of the loops bounds a disc, then this routine might (but is not
        guaranteed to) raise BoundsDisc.

        This is equivalent to calling loop.shorten() for each loop in self.

        Pre-condition:
        --> Each loop in self must be an instance of a subclass of
            EmbeddedLoop that implements the shorten() routine.

        Returns:
            True if and only if the union of embedded loops was successfully
            shortened.
        """
        changed = False
        for embLoop in self:
            if embLoop.shorten():
                changed = True
        return

    def _minimiseBoundaryImpl(self):
        """
        Ensures that the ambient triangulation has the smallest possible
        number of boundary triangles, potentially adding tetrahedra to do
        this.

        A side-effect of calling this routine is that it will shorten the
        union of embedded loops if possible.

        The default implementation of this routine relies on the following
        subroutines, not all of which are implemented by default:
        --> shorten()
        --> _findBoundaryMove()
        Thus, subclasses that require this routine must either:
        --> override this routine; or
        --> ensure that all of the aforementioned subroutines are suitably
            implemented.

        In the default implementation, every call to _findBoundaryMove() is
        immediately preceded by a call to shorten(). Thus, for subclasses that
        use the default implementation of this routine, any post-conditions of
        shorten() can be assumed as pre-conditions of _findBoundaryMove().

        If one of the embedded loops bounds a disc, then this routine might
        (but is not guaranteed to) raise BoundsDisc.

        The following are guaranteed to hold once this routine is finished:
        --> Every 2-sphere boundary component will have exactly two triangles
            and three vertices.
        --> Every projective plane boundary component will have exactly two
            triangles and two vertices.
        --> Every other boundary component will have exactly one vertex.

        The changes that this routine performs can always be expressed using
        only the following operations:
        --> Shortening some loop by redirecting it across triangular faces.
        --> Close book moves and/or layerings on self.triangulation().
        In particular, this routine never creates new vertices, and it never
        creates a non-vertex-linking normal disc or 2-sphere if there was not
        one before.

        If the ambient triangulation is currently oriented, then this routine
        guarantees to preserve the orientation.

        Adapted from Regina's Triangulation3.minimiseBoundary().

        Precondition:
        --> The ambient triangulation (i.e., self.triangulation()) is valid.

        Returns:
            True if and only if this triangulation with embedded loops was
            changed. In other words, a return value of False indicates that:
            (1) the union of loops could not be shortened; and
            (2) every boundary component of the ambient triangulation was
                already minimal to begin with.
        """
        changed = False
        while True:
            # Shorten the loops to minimise the number of special cases.
            if self.shorten():  # Might raise BoundsDisc.
                changed = True

            # Is there a move we can perform to reduce the number of boundary
            # triangles?
            moveDetails = self._findBoundaryMove()  # Might raise BoundsDisc.
            if moveDetails is None:
                # At this point, all boundary components must be minimal.
                return changed

            # Perform the move, and then update this loop.
            #
            # The _findBoundaryMove() routine should already ensure legality
            # of the close book move, so no need to check before performing.
            changed = True
            edge, doLayer, edgeEmbeddingsData = moveDetails
            if doLayer:
                edge = layerOn(edge).edge(5)
            self._tri.closeBook( edge, False, True )
            self.setFromEdgeEmbeddings(edgeEmbeddingsData)
        return

    def _findBoundaryMove(self):
        """
        Returns details of a boundary move that simplifies the boundary of
        self.triangulation(), or None if the boundary is already minimal.

        In detail, in the case where the boundary is not yet minimal, this
        routine guarantees to find a move that reduces the number of boundary
        triangles by two (without changing the topology of the union of
        embedded loops). The return value will be a tuple that describes this
        move using the following data:
        (0) A boundary edge e on which to perform the move.
        (1) A boolean indicating whether we need to layer across e. If this
            is True, then the move we perform will be to first layer across
            e, and then perform a close book move on the newly layered edge.
            Otherwise, the move will simply be a close book move on e.
        (2) Data for reconstructing the union of embedded loops in the new
            triangulation that results from this move, as required by the
            setFromEdgeEmbeddings() routine.

        The TriangulationWithEmbeddedLoops base class does not implement this
        routine, so subclasses that require this routine must provide an
        implementation.

        If necessary, an implementation may raise BoundsDisc when it detects
        an embedded loop that bounds a disc.

        Also, an implementation is allowed to assume, as a pre-condition, that
        all embedded loops have already been shortened as much as possible.
        """
        raise NotImplementedError()

    #TODO WORKING HERE

    #TODO Document precisely how strongly we can guarantee minimality.
    def _minimiseVerticesImpl(self):
        """
        Ensures that this triangulation with embedded loops has the smallest
        possible number of vertices, potentially adding tetrahedra to do this.

        The default implementation of this routine relies on the following
        subroutines, not all of which are implemented by default:
        --> shorten()
        --> _minimiseBoundaryImpl()
        --> _findSnapEdge()
        Thus, subclasses that require this routine must either:
        --> override this routine; or
        --> ensure that all of the aforementioned subroutines are suitably
            implemented.

        A side-effect of calling this routine is that it will shorten the
        union of embedded loops if possible.

        If one of the embedded loops bounds a disc, then this routine might
        (but is not guaranteed to) raise BoundsDisc.

        Subclasses may decide precisely what it means for the number of
        vertices to be minimised. Here are some typical examples of conditions
        that might constitute minimality:
        --> If the ambient triangulation is closed, then it has precisely one
            vertex per embedded loop.
        --> If the ambient triangulation has real boundary, and every embedded
            loop is either (entirely) internal or (entirely) boundary, then:
            --- there is exactly one internal vertex for each internal loop;
            --- no boundary loop is incident to a 2-sphere or projective plane
                boundary component;
            --- every 2-sphere boundary component has exactly two triangles
                and three vertices;
            --- every projective plane boundary component has exactly two
                triangles and two vertices; and
            --- the total number of vertices in all the other boundary
                components is equal to the number of boundary loops.

        The changes that this routine performs can always be expressed using
        only the following operations:
        --> Shortening some loop by redirecting it across triangular faces.
        --> Close book moves, layerings, and/or snap edge moves on
            self.triangulation().
        In particular, this routine never creates new vertices.

        If the ambient triangulation is currently oriented, then this routine
        guarantees to preserve the orientation.

        Adapted from Regina's Triangulation3.minimiseVertices().

        Precondition:
        --> The ambient triangulation (i.e., self.triangulation()) is valid.

        Returns:
            True if and only if this triangulation with embedded loops was
            changed. In other words, a return value of False indicates that:
            (1) the union of loops could not be shortened; and
            (2) the number of vertices in the ambient triangulation was
                already minimal to begin with.
        """
        # Start by minimising the boundary.
        changed = self._minimiseBoundaryImpl()  # Might raise BoundsDisc.

        # All that remains now is to remove unnecessary internal vertices.
        # We do not currently have an implementation of collapseEdge() that
        # keeps track of how edges get relabelled after the move, so we rely
        # entirely on the snapEdge() routine.
        while True:
            # Shorten the union of loops to minimise the number of special
            # cases.
            if self.shorten():  # Might raise BoundsDisc.
                changed = True

            # Is there a snap edge move we can perform to reduce the number
            # of vertices?
            moveDetails = self._findSnapEdge()
            if moveDetails is None:
                # At this point, there are no more unnecessary internal
                # vertices.
                return changed
            changed = True
            edge, edgeEmbeddingsData = moveDetails

            # Perform the snap, and then update this ideal loop. The
            # _findSnapEdge() routine should ensure that the snap is legal, so
            # no need to check this again.
            snapEdge( edge, False, True )
            self.setFromEdgeEmbeddings(edgeEmbeddingsData)
        return

    #TODO Check pre-conditions and exceptions.
    def _findSnapEdge(self):
        """
        Returns details of a snap edge move that can be used to reduce the
        number of vertices in self.triangulation(), or None if the number of
        vertices is already minimal.

        As in the documentation for _minimiseVerticesImpl(), subclasses may
        decide on precisely what minimal should mean.

        In detail, in the case where the number of vertices is not yet
        minimal, this routine returns a tuple consisting of the following:
        (0) An edge on which a snap edge move can be performed.
        (1) Data for reconstructing the union of embedded loops in the new
            triangulation that results from this move, as required by the
            setFromEdgeEmbeddings() routine.

        The TriangulationWithEmbeddedLoops base class does not implement this
        routine, so subclasses that require this routine must provide an
        implementation.

        Pre-condition:
        --> The union of loops cannot be shortened.
        --> If the ambient triangulation has real boundary, then this
            boundary has already been minimised.
        """
        raise NotImplementedError()

    #TODO
    pass


class EdgeIdealTriangulation(TriangulationWithEmbeddedLoops):
    """
    An edge-ideal triangulation; that is, a 3-manifold triangulation
    containing a disjoint union of IdealLoop objects in its interior.

    Some of the routines provided by this class might fail if one of the ideal
    loops bounds an embedded disc in the ambient triangulation (though these
    routines might nevertheless succeed in spite of the existence of such a
    disc). This class raises BoundsDisc whenever such a failure occurs.

    This class has two core features:
    (1) It provides methods to simplify the ambient 3-manifold triangulation,
        while ensuring that the topological embedding of the union of ideal
        loops is always preserved.
    (2) It acts as a container of IdealLoop objects, which are indexed in an
        arbitrary order (but the order is kept consistent no matter how much
        the ambient triangulation is simplified). In detail, for any instance
        triLoops of this class:
        --> (el in triLoops) is True if and only if el is one of the IdealLoop
            objects that is embedded in the ambient triangulation (note that
            equality of IdealLoop objects is determined by their location in
            memory, and so for instance clones will not be considered equal).
        --> len(triLoops) is the number of IdealLoop objects that are embedded
            in the ambient triangulation
        --> iterating through triLoops yields all the IdealLoop objects that
            are embedded in the ambient triangulation, in the order in which
            they are indexed
        --> for i between 0 and (len(triLoops) - 1), inclusive, triLoops[i]
            returns the IdealLoop at index i
    """
    # Base class for loops embedded in the ambient triangulation.
    _LOOP_CLASS = IdealLoop

    def __init__( self, loops ):
        """
        Creates an edge-ideal triangulation containing the disjoint union of
        the given collection of ideal loops.

        Precondition:
        --> loops is nonempty.
        --> The elements of loops are all IdealLoop objects lying disjointly
            in the interior of the same ambient 3-manifold triangulation.
        """
        super().__init__(loops)
        return

    def shorten(self):
        """
        Shortens the union of ideal loops.

        In detail, if some ideal loop meets any triangle F in exactly two
        distinct edges, then it can be shortened by replacing these two edges
        with the third edge of F.

        This routine performs such shortenings until no further shortening is
        possible. If at least one such shortening occurred, then this routine
        will return True. Otherwise, this routine will leave the union of
        ideal loops entirely untouched, and will return False.

        This routine raises BoundsDisc if self.triangulation() includes a
        triangular face F that forms an embedded disc whose boundary is given
        by one of the ideal loops.

        Returns:
            True if and only if the union of ideal loops was successfully
            shortened.
        """
        # IdealLoop provides an appropriate implementation of shorten(), so we
        # can just use the default implementation from the base class.
        return super().shorten()

    def minimiseBoundary(self):
        """
        Ensures that the ambient triangulation has the smallest possible
        number of boundary triangles, potentially adding tetrahedra to do
        this.

        A side-effect of calling this routine is that it will shorten the
        ideal loops if possible.

        If some ideal loop bounds a disc, then this routine might (but is not
        guaranteed to) raise BoundsDisc.

        The following are guaranteed to hold once this routine is finished:
        --> Every 2-sphere boundary component will have exactly two triangles
            and three vertices.
        --> Every projective plane boundary component will have exactly two
            triangles and two vertices.
        --> Every other boundary component will have exactly one vertex.

        The changes that this routine performs can always be expressed using
        only the following operations:
        --> Shortening a loop by redirecting it across triangular faces.
        --> Close book moves and/or layerings on self.triangulation().
        In particular, this routine never creates new vertices, and it never
        creates a non-vertex-linking normal disc or 2-sphere if there was not
        one before.

        If the ambient triangulation is currently oriented, then this routine
        guarantees to preserve the orientation.

        Adapted from Regina's Triangulation3.minimiseBoundary().

        Precondition:
        --> The ambient triangulation (i.e., self.triangulation()) is valid.

        Returns:
            True if and only if this edge-ideal triangulation was changed. In
            other words, a return value of False indicates that:
            (1) the union of loops could not be shortened; and
            (2) every boundary component of the ambient triangulation was
                already minimal to begin with.
        """
        # Can use the default implementation provided we supply an
        # implementation for _findBoundaryMove().
        return self._minimiseBoundaryImpl()

    def _findBoundaryMove(self):
        # Precondition:
        #   --> The union of loops cannot be shortened.

        # Find a boundary component that is not yet minimal.
        for bc in self._tri.boundaryComponents():
            if bc.countTriangles() <= 2 or bc.countVertices() <= 1:
                # Current boundary component is already minimal.
                continue

            # First try to find a close book move, which does not increase
            # the number of tetrahedra.
            for edge in bc.edges():
                if self._tri.closeBook( edge, True, False ):
                    return ( edge,
                            False,  # Close book w/out layering.
                            self._edgeEmbeddingsData() )

            # We could not find a close book move.
            # In this case, because bc is non-minimal, there must be a
            # boundary edge e that joins two distinct vertices, and we can
            # simplify bc by layering across e and then performing a close
            # book move on the newly layered edge.
            for edge in bc.edges():
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
                        True,   # Layer before doing close book.
                        self._edgeEmbeddingsData() )

            # We should never reach this point.
            raise AssertionError(
                    "_findBoundaryMove() failed unexpectedly." )

        # If we fell out of the boundary component loop, then all boundary
        # components are minimal.
        return None

    def minimiseVertices(self):
        """
        Ensures that this edge-ideal triangulation has the smallest possible
        number of vertices, potentially adding tetrahedra to do this.

        A side-effect of calling this routine is that it will shorten the
        ideal loops if possible.

        If some ideal loop bounds a disc, then this routine might (but is not
        guaranteed to) raise BoundsDisc.

        The following are guaranteed to hold once this routine is finished:
        --> If the ambient triangulation is closed, then it will have
            precisely one (internal) vertex per ideal loop.
        --> If the ambient triangulation has real boundary, then:
            --- there will be exactly one internal vertex for each ideal loop;
            --- every 2-sphere boundary component will have exactly two
                triangles and three vertices;
            --- every projective plane boundary component will have exactly
                two triangles and two vertices; and
            --- every other boundary component will have exactly one vertex.

        The changes that this routine performs can always be expressed using
        only the following operations:
        --> Shortening some loop by redirecting it across triangular faces.
        --> Close book moves, layerings, and/or snap edge moves on
            self.triangulation().
        In particular, this routine never creates new vertices.

        If the ambient triangulation is currently oriented, then this routine
        guarantees to preserve the orientation.

        Adapted from Regina's Triangulation3.minimiseVertices().

        Precondition:
        --> The ambient triangulation (i.e., self.triangulation()) is valid.

        Returns:
            True if and only if this edge-ideal triangulation was changed. In
            other words, a return value of False indicates that:
            (1) the union of loops could not be shortened; and
            (2) the number of vertices in the ambient triangulation was
                already minimal to begin with.
        """
        # Can use the default implementation provided we supply an
        # implementation for _findSnapEdge().
        return self._minimiseVerticesImpl()

    def _findSnapEdge(self):
        # Pre-condition:
        #   --> The union of loops cannot be shortened.
        #   --> The boundary of self._tri has been minimised.

        # Find a suitable edge on which to perform a snap edge move. We
        # minimise the number of special cases by prioritising edges that
        # belong to one of the ideal loops.
        for iloop in self:
            if len(iloop) == 1:
                continue

            # We can shorten iloop by snapping of its edges. Here, we
            # choose the last edge of iloop.
            data = self._edgeEmbeddingsData( remove={ iloop[-1] } )
            return ( self._tri.edge( iloop[-1] ), data )

        # At this point, the boundary has been minimised, and every ideal loop
        # already has length 1. However, there might still be some other
        # internal vertices that we can remove.
        for edge in self._tri.edges():
            if not snapEdge( edge, True, False ):
                # Snap edge is not legal.
                continue

            # The snap edge move is legal, but we only want to perform this
            # move if it will remove an internal vertex that is not incident
            # to any of the ideal loops.
            for i in range(2):
                v = edge.vertex(i)
                if ( v.isBoundary() or
                    ( v.index() in self.loopVertexIndices() ) ):
                    continue

                # This edge has an endpoint which is both internal and not
                # incident to any ideal loops. Snapping this edge will
                # effectively merge this vertex with whatever vertex is at the
                # opposite endpoint.
                return ( edge, self._edgeEmbeddingsData() )

        return

    #TODO
    pass


#TODO We never use 2-sphere or projective plane boundary components, and
#       handling them is a mess in terms of special cases (especially for the
#       minimiseBoundary() and minimiseVertices() routines). Should we just
#       exclude such boundary components?
class TriangulationWithBoundaryLoops(TriangulationWithEmbeddedLoops):
    """
    A 3-manifold triangulation containing up to one BoundaryLoop per boundary
    component.

    Some of the routines provided by this class might fail if one of the
    boundary loops bounds an embedded disc in the ambient triangulation
    (though these routines might nevertheless succeed in spite of the
    existence of such a disc). This class raises BoundsDisc whenever such a
    failure occurs.

    This class has two core features:
    (1) It provides methods to simplify the ambient 3-manifold triangulation,
        while ensuring that the topological embedding of the union of boundary
        loops is always preserved.
    (2) It acts as a container of BoundaryLoop objects, which are indexed in
        an arbitrary order (but the order is kept consistent no matter how
        much the ambient triangulation is simplified). In detail, for any
        instance triLoops of this class:
        --> (el in triLoops) is True if and only if el is one of the
            BoundaryLoop objects that is embedded in the ambient triangulation
            (note that equality of BoundaryLoop objects is determined by their
            location in memory, and so for instance clones will not be
            considered equal).
        --> len(triLoops) is the number of BoundaryLoop objects that are
            embedded in the ambient triangulation
        --> iterating through triLoops yields all the BoundaryLoop objects
            that are embedded in the ambient triangulation, in the order in
            which they are indexed
        --> for i between 0 and (len(triLoops) - 1), inclusive, triLoops[i]
            returns the BoundaryLoop at index i
    """
    # Base class for loops embedded in the ambient triangulation.
    _LOOP_CLASS = BoundaryLoop

    def __init__( self, loops ):
        """
        Creates a triangulation containing the given union of boundary loops.

        Precondition:
        --> loops is nonempty.
        --> The elements of loops are all BoundaryLoop objects lying on
            distinct boundary components of the same ambient 3-manifold
            triangulation.
        """
        super().__init__(loops)
        return

    def shorten(self):
        """
        Shortens the union of boundary loops.

        In detail, if some boundary loop meets any boundary triangle F in
        exactly two distinct edges, then it can be shortened by replacing
        these two edges with the third edge of F.

        This routine performs such shortenings until no further shortening is
        possible. If at least one such shortening occurred, then this routine
        will return True. Otherwise, this routine will leave this
        triangulation with boundary loops entirely untouched, and will return
        False.

        This routine raises BoundsDisc if self.triangulation() includes a
        boundary triangular face F that forms an embedded disc whose boundary
        is given by one of the boundary loops.

        Returns:
            True if and only if the union of boundary loops was successfully
            shortened.
        """
        # BoundaryLoop provides an appropriate implementation of shorten(), so
        # we can just use the default implementation from the base class.
        return self.shorten()

    def minimiseBoundary(self):
        """
        Ensures that the ambient triangulation has the smallest possible
        number of boundary triangles, potentially adding tetrahedra to do
        this.

        A side-effect of calling this routine is that it will shorten the
        boundary loops if possible.

        If some boundary loop bounds a disc, then this routine might (but is
        not guaranteed to) raise BoundsDisc.

        The following are guaranteed to hold once this routine is finished:
        --> Every 2-sphere boundary component will have exactly two triangles
            and three vertices.
        --> Every projective plane boundary component will have exactly two
            triangles and two vertices.
        --> Every other boundary component will have exactly one vertex.

        The changes that this routine performs can always be expressed using
        only the following operations:
        --> Shortening a loop by redirecting it across triangular boundary
            faces.
        --> Close book moves and/or layerings on self.triangulation().
        In particular, this routine never creates new vertices, and it never
        creates a non-vertex-linking normal disc or 2-sphere if there was not
        one before.

        If the ambient triangulation is currently oriented, then this routine
        guarantees to preserve the orientation.

        If some boundary loop bounds a disc, then this routine might (but is
        not guaranteed to) raise BoundsDisc.

        Adapted from Regina's Triangulation3.minimiseBoundary().

        Precondition:
        --> The ambient triangulation (i.e., self.triangulation()) is valid.

        Returns:
            True if and only if this triangulation with boundary loops was
            changed. In other words, a return value of False indicates that:
            (1) the union of loops could not be shortened; and
            (2) every boundary component of the ambient triangulation was
                already minimal to begin with.
        """
        # Can use the default implementation provided we supply an
        # implementation for _findBoundaryMove().
        return self._minimiseBoundaryImpl()

    #TODO Is it worthwhile to be more careful with the implementation, and
    #       thereby allow multiple boundary loops per component?
    def _findBoundaryMove(self):
        # Exceptions:
        #   --> Might raise BoundsDisc.
        #
        # Precondition:
        #   --> The union of loops cannot be shortened.

        # Prioritise moves that reduce the length of the boundary loops. If
        # possible, use close book moves so that we do not introduce too many
        # new tetrahedra.
        for bloop in self:
            if len(bloop) == 1:
                continue
            if bloop.boundaryComponent().countTriangles() == 2:
                # We have a minimal boundary component bc with more than one
                # vertex (because bloop already has more than one vertex).
                # Thus bc must be either a 2-sphere or a projective plane.
                # Because bloop cannot be shortened, it is incident to each
                # triangle of bc in either one edge or three edges.
            #TODO
            if ( len(bloop) == 1 or
                bloop.boundaryComponent().countTriangles() == 2 ):
                continue

            # We have a boundary loop in a non-minimal boundary component.
            # First try to find a close book move.
            if len(bloop) > 2:
                for edge in bloop.boundaryComponent().edges():
                    # Check eligibility of close book move, but do *not*
                    # perform yet.
                    if not bloop.triangulation.closeBook(
                            edge, True, False ):
                        continue

                    # Because we have assumed that bloop cannot be
                    # shortened, the only way this close book move can
                    # reduce the length of bloop is if bloop meets either:
                    #   --> both left edges in the diagram below, or
                    #   --> both right edges.
                    #
                    #                   2
                    #                   •
                    #                  / \ front
                    #                 /   \
                    #               0•-----•1
                    #               0•-----•1
                    #                 \   /
                    #                  \ / back
                    #                   •
                    #                   3
                    #
                    ftet = edge.front().tetrahedron()
                    fver = edge.front().vertices()
                    btet = edge.back().tetrahedron()
                    bver = edge.back().vertices()
                    for v in range(2):
                        fei = ftet.edge( fver[v], fver[2] ).index()
                        bei = btet.edge( bver[v], bver[3] ).index()
                        if ( fei in bloop ) and ( bei in bloop ):
                            # After performing the close book move, the
                            # loop will no longer include edges fei and
                            # bei.
                            data = self._edgeEmbeddingsData(
                                    remove={ fei, bei } )
                            return ( edge,
                                    False,  # Close book w/out layering.
                                    data )

            # Resort to layering a new tetrahedron to facilitate a close
            # book move that effectively removes one of the edges of
            # bloop. This operation is guaranteed to be legal for any edge
            # of bloop; here, we choose to perform it on the last edge.
            data = self._edgeEmbeddingsData( remove={ bloop[-1] } )
            return ( self._tri.edge( bloop[-1] ),
                    True,   # Layer before doing close book.
                    data )

        # At this point, for each boundary loop bloop, it must already be the
        # case that either:
        #   --> bloop has length 1; or
        #   --> if this is impossible, which can only be the case if bloop
        #       lies in a 2-sphere or projective plane boundary component,
        #       then we at least have that bloop lies inside a two-triangle
        #       (hence minimal) boundary component.
        # Our goal is to minimise any remaining non-minimal boundary
        # components without touching the (necessarily length-1) loops that
        # they contain.
        avoid = self.loopEdgeIndices()
        for bc in self._tri.boundaryComponents():
            if bc.countTriangles() == 2 or bc.countVertices() == 1:
                # Already minimal.
                continue

            # First try to find a close book move.
            for edge in bc.edges():
                if edge.index() in avoid:
                    continue
                if self._tri.closeBook( edge, True, False ):
                    return ( edge,
                            False,  # Close book w/out layering.
                            self._edgeEmbeddingsData() )

            # We could not find a suitable close book move, so our plan now is
            # to find a boundary edge e that joins two distinct vertices (such
            # an edge e must exist). We can layer over such an edge e, and
            # then do a close book move on the newly-layered edge.
            for edge in bc.edges():
                if edge.vertex(0) == edge.vertex(1):
                    # Note that this automatically avoids any edges involved
                    # in the boundary loops, since we know that bc contains at
                    # most a single such loop, and this loop will have length
                    # 1 if it exists.
                    continue

                # The layering is illegal if this edge is incident to the same
                # boundary triangle F on both sides (rather than two distinct
                # triangles). But in that scenario, F forms a disc with
                # boundary given by the third edge b of F. There must be a
                # close book move available on b, so if we didn't perform this
                # move already then it must be the case that b forms the loop
                # embedded in this boundary component.
                front = edge.front()
                back = edge.back()
                if ( front.tetrahedron().triangle( front.vertices()[3] ) ==
                        back.tetrahedron().triangle( back.vertices()[2] ) ):
                    raise BoundsDisc()
                else:
                    return ( edge,
                            True,   # Layer before performing close book.
                            self._edgeEmbeddingsData() )

            # We should never reach this point.
            raise AssertionError(
                    "_findBoundaryMove() failed unexpectedly." )

        # If we fell out of the boundary component loop, then all boundary
        # components are minimal.
        return None

    def minimiseVertices(self):
        """
        Ensures that this triangulation with boundary loops has the smallest
        possible number of vertices, potentially adding tetrahedra to do this.

        A side-effect of calling this routine is that it will shorten the
        boundary loops if possible.

        If some boundary loop bounds a disc, then this routine might (but is
        not guaranteed to) raise BoundsDisc.

        The following are guaranteed to hold once this routine is finished:
        --> There will be no internal vertices.
        --> Every 2-sphere boundary component will have exactly two triangles
            and three vertices.
        --> Every projective plane boundary component will have exactly two
            triangles and two vertices.
        --> Every other boundary component will have exactly one vertex.

        The changes that this routine performs can always be expressed using
        only the following operations:
        --> Shortening some loop by redirecting it across triangular faces.
        --> Close book moves, layerings, and/or snap edge moves on
            self.triangulation().
        In particular, this routine never creates new vertices.

        If the ambient triangulation is currently oriented, then this routine
        guarantees to preserve the orientation.

        Adapted from Regina's Triangulation3.minimiseVertices().

        Precondition:
        --> The ambient triangulation (i.e., self.triangulation()) is valid.

        Returns:
            True if and only if this triangulation with boundary loops was
            changed. In other words, a return value of False indicates that:
            (1) the union of loops could not be shortened; and
            (2) the number of vertices in the ambient triangulation was
                already minimal to begin with.
        """
        # Can use the default implementation provided we supply an
        # implementation for _findSnapEdge().
        return self._minimiseVerticesImpl()

    def _findSnapEdge(self):
        # Pre-condition:
        #   --> The union of loops cannot be shortened.
        #   --> The boundary of self._tri has been minimised.

        # Find a suitable edge on which to perform a snap edge move (just
        # check whether the move is legal, do not perform yet).
        for edge in self._tri.edges():
            if snapEdge( edge, True, False ):
                # The loops should all lie entirely in the boundary, so the
                # snap edge move should never change the loops.
                return ( edge, self._edgeEmbeddingsData() )
        return

    #TODO
    pass


#TODO Test suite.
