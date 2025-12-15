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
#TODO Reimplement all the simplification methods so that they can handle
#   unions of more than one embedded loop.


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

    def _shortenImpl(self):
        """
        Shortens the union of embedded loops by looking for triangles that
        intersect a loop in two edges, and redirecting this loop to use the
        third edge.

        If no shortening is possible, then this triangulation with embedded
        loops will remain entirely untouched.

        The default implementation of this routine requires the helper routine
        self._redirectCandidates(), which has the following pre-condition:
        --> For every loop in self, loop must be a subclass of EmbeddedLoop
            that implements the loop._redirectCandidates() routine.

        If one of the loops bounds a disc, then this routine might (but is not
        guaranteed to) raise BoundsDisc.

        Returns:
            True if and only if the union of embedded loops was successfully
            shortened.
        """
        if self.countLoopEdges() == len(self):
            # Every loop already uses just a single edge.
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
        _shortenImpl() routine should attempt to redirect the union of
        embedded loops.

        The default implementation of this routine has the following
        pre-condition:
        --> For every loop in self, loop must be a subclass of EmbeddedLoop
            that implements the loop._redirectCandidates() routine.
        """
        for embLoop in self:
            for candidate in embLoop._redirectCandidates():
                yield candidate
        return

    def _attemptRedirect( self, face ):
        r"""
        Attempts to redirect one of the embedded loops across the given face.

        If the given face intersects a loop in exactly two *distinct* edges,
        then we can redirect this loop along the third edge of the given face.

                    Before redirect         After redirect
                           •                       •
                          / \
                         /   \
                        /     \
                       •       •               •-------•

        If such a redirect is possible, then this routine performs the
        redirect and returns True. Otherwise, this routine leaves this
        triangulation with embedded loops entirely untouched and returns
        False.

        Parameters
        --> face    the triangular face across which to attempt a redirect of
                    some loop

        Returns:
            True if and only if this routine successfully performs the
            requested redirect
        """
        # Because all the loops are disjointly embedded, we can independently
        # run each individual loop's own implementation of _attemptRedirect().
        for embLoop in self:
            if embLoop._attemptRedirect(face):
                return True
        return False

    #TODO WORKING HERE

    #TODO Make sure that the documentation matches the new implementation.
    def _minimiseBoundaryImpl(self):
        """
        Ensures that the ambient triangulation has the smallest possible
        number of boundary triangles, potentially adding tetrahedra to do
        this.

        The default implementation of this routine requires the following
        helper routines:
        --> _shortenImpl()
        --> _findBoundaryMove()
        Thus, subclasses that require this routine must either:
        --> override this routine; or
        --> ensure that the aforementioned helper routines are suitably
            implemented.

        A side-effect of calling this routine is that it will shorten the
        union of embedded loops if possible.

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
            if self._shortenImpl():     # Might raise BoundsDisc.
                changed = True

            # Is there a move we can perform to reduce the number of boundary
            # triangles?
            moveDetails = self._findBoundaryMove()
            if moveDetails is None:
                # At this point, all boundary components must be minimal.
                return changed

            # Perform the move, and then update this loop.
            #
            # The _findBoundaryMove() routine should already ensure legality
            # of the close book move, so no need to check before performing.
            changed = True
            edge, doLayer, newLoopData = moveDetails
            if doLayer:
                edge = layerOn(edge).edge(5)
            self._tri.closeBook( edge, False, True )
            self.setFromEdgeEmbeddings(newLoopData)
        return

    #TODO Document exceptions and pre-conditions.
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

        If some ideal loop bounds a disc, then this routine might (but is not
        guaranteed to) raise BoundsDisc.

        Returns:
            True if and only if the union of ideal loops was successfully
            shortened.
        """
        # IdealLoop provides an appropriate implementation of
        # _redirectCandidates(), so we can just use the default implementation
        # of shortening
        return self._shortenImpl()  # Might raise BoundsDisc.

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
            raise RuntimeError(
                    "_findBoundaryMove() failed unexpectedly." )

        # If we fell out of the boundary component loop, then all boundary
        # components are minimal.
        return None

    #TODO
    pass


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

        If some boundary loop bounds a disc, then this routine might (but is
        guaranteed to) raise BoundsDisc.

        Returns:
            True if and only if the union of boundary loops was successfully
            shortened.
        """
        # BoundaryLoop provides an appropriate implementation of
        # _redirectCandidates(), so we can just use the default implementation
        # of shortening
        return self._shortenImpl()  # Might raise BoundsDisc.

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

    #TODO Update to:
    #       --> handle multiple loops, and
    #       --> track orientations.
    def _findBoundaryMove(self):
        # Exceptions:
        #   --> Might raise BoundsDisc.
        #
        # Precondition:
        #   --> The union of loops cannot be shortened.

        # Prioritise moves that reduce the length of the boundary loops. If
        # possible, use close book moves so that we do not introduce too many
        # new tetrahedra.
        if self.countLoopEdges() > len(self):
            for bloop in self:
                if ( len(bloop) == 1 or
                    bloop.boundaryComponent().countTriangles() == 2 ):
                    continue

                #TODO Maybe outsource some of this to BoundaryLoop.

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

        #TODO Re-implement.

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

    #TODO
    pass


#TODO Test suite.
