"""
3-manifold triangulations containing disjoint unions embedded loops.

The classes here provide methods for simplifying the ambient 3-manifold
triangulation, while preserving the topological embedding of the loops.
"""
from regina import *
from loop import EmbeddedLoop, IdealLoop, BoundaryLoop
from aux.looperror import BoundsDisc
from aux.edgeemb import embeddingsFromEdgeIndices
from aux.surface import SurfaceType, hasOnlyNonTrivialBoundaryCurves
from retriangulate.moves import twoThree, threeTwo, twoZero, twoOne, fourFour
from retriangulate.insert import snapEdge, layerOn
from retriangulate.edgelabel import EdgeLabelling


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

        The given loops can be any iterable that yields EmbeddedLoop objects.
        For example, it could be a Python list of EmbeddedLoop objects, or it
        could be a TriangulationWithEmbeddedLoops object.

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

    def _edgeLab(self):
        """
        Returns an EdgeLabelling that only tracks the edges involved in the
        embedded loops.
        """
        return EdgeLabelling(
                self._tri,
                { ei: self._tri.edge(ei).front()
                 for ei in self.loopEdgeIndices() } )

    def _setFromRelab( self, relab ):
        """
        Sets this triangulation with embedded loops using the relabelling
        described by the given EdgeLabelling.

        This routine is for internal use only. The purpose of this routine is
        to update the embedded loops whenever the ambient triangulation has
        been modified by a local move. See the twoThree, threeTwo, twoZero,
        twoOne, and fourFour routines from retriangulate/moves.py for examples
        of how relabellings are specified.

        Pre-condition:
        --> The given EdgeLabelling relab tracks every index ei in
            self.loopEdgeIndices().
        """
        for embLoop in self:
            embLoop._setFromRelab(relab)
        return

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

    def isLoopEdgeIndex( self, edgeIndex ):
        """
        Does the given edge index correspond to an edge that belongs to one of
        the embedded loops?
        """
        for embLoop in self:
            if edgeIndex in embLoop:
                return True
        return False

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

    #TODO Which of the following do we really need?
    #       --> intersects()
    #       --> incidentLoopWeights()
    #       --> weight()
    #       --> incidentLoopIndices()

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

    def incidentLoopWeights( self, surf ):
        """
        Returns a dictionary mapping loop indices to their weights with
        respect to the given normal surface.

        Only loops with positive weights will be included in the returned
        dictionary.

        Precondition:
        --> The given normal surface is embedded in self.triangulation().
        """
        ans = dict()
        for i, embLoop in enumerate(self):
            wt = embLoop.weight(surf)
            if wt > 0:
                ans[i] = wt
        return ans

    def weight( self, surf ):
        """
        Returns the number of times the union of embedded loops intersects
        the given normal surface surf.

        Precondition:
        --> The given normal surface is embedded in self.triangulation().
        """
        return sum( self.incidentLoopWeights(surf).values() )

    def incidentLoopIndices( self, surf ):
        """
        Returns a set consisting of the indices of the embedded loops that
        are incident to the given normal surface.

        Precondition:
        --> The given normal surface is embedded in self.triangulation().
        """
        return set( self.incidentLoopWeights(surf).keys() )

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

    def minimiseBoundary(self):
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

        Subclasses may decide precisely what it means for the number of
        boundary triangles to be minimised. Typically, the boundary will be
        minimised when the following conditions hold:
        --> Every 2-sphere boundary component has exactly two triangles and
            three vertices.
        --> Every projective plane boundary component has exactly two
            triangles and two vertices.
        --> Every other boundary component has exactly one vertex.

        For the default implementation, the minimum number of boundary
        triangles is achieved precisely when the _findBoundaryMove() routine
        returns None.

        If no exceptions are raised, then the following are guaranteed to hold
        once this routine terminates:

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
            # The _findBoundaryMove() routine ensures legality of the close
            # book move.
            changed = True
            edge, doLayer, edgeEmbeddingsData = moveDetails
            if doLayer:
                edge = layerOn(edge).edge(5)
            self._tri.closeBook(edge)
            self.setFromEdgeEmbeddings(edgeEmbeddingsData)
        return

    def _findBoundaryMove(self):
        """
        Returns details of a boundary move that simplifies the boundary of
        self.triangulation(), or None if the boundary is already minimal.

        As in the documentation for minimiseVertices(), subclasses may decide
        on precisely what minimal should mean.

        In the case where the boundary is not yet minimal, this routine
        guarantees to find a move that reduces the number of boundary
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

        Also, an implementation is allowed to assume any post-conditions of
        the shorten() routine as pre-conditions for this routine.
        """
        raise NotImplementedError()

    def minimiseVertices(self):
        """
        Ensures that this triangulation with embedded loops has the smallest
        possible number of vertices, potentially adding tetrahedra to do this.

        A side-effect of calling this routine is that it will shorten the
        union of embedded loops if possible.

        The default implementation of this routine relies on the following
        subroutines, not all of which are implemented by default:
        --> shorten()
        --> minimiseBoundary()
        --> _findSnapEdge()
        Thus, subclasses that require this routine must either:
        --> override this routine; or
        --> ensure that all of the aforementioned subroutines are suitably
            implemented.

        In the default implementation, every call to _findSnapEdge() is
        immediately preceded by a call to shorten(). The default
        implementation also begins by calling minimiseBoundary(), and then
        never modifies the boundary again. Thus, for subclasses that use the
        default implementation of this routine, any post-conditions of the
        following subroutines can be assumed as pre-conditions of
        _findSnapEdge():
        --> shorten()
        --> minimiseBoundary()

        If one of the embedded loops bounds a disc, then this routine might
        (but is not guaranteed to) raise BoundsDisc.

        Subclasses may decide precisely what it means for the number of
        vertices to be minimised. For example, if the ambient triangulation is
        closed, then the minimum possible number of vertices would be equal to
        the number of embedded loops.

        For the default implementation, the minimum number of vertices is
        achieved precisely when the _findSnapEdge() routine returns None.

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
        changed = self.minimiseBoundary()   # Might raise BoundsDisc.

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

    def _findSnapEdge(self):
        """
        Returns details of a snap edge move that can be used to reduce the
        number of vertices in self.triangulation(), or None if the number of
        vertices is already minimal.

        As in the documentation for minimiseVertices(), subclasses may decide
        on precisely what minimal should mean.

        In the case where the number of vertices is not yet minimal, this
        routine returns a tuple consisting of the following:
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

    def _simplifyMonotonicImpl( self, include32 ):
        """
        Uses 2-0 edge, 2-1 edge, and (optionally) 3-2 moves to monotonically
        reduce the number of tetrahedra in the ambient triangulation, while
        leaving the embedded loops untouched.

        If one of the embedded loops bounds a disc, then this routine might
        (but is not guaranteed to) raise BoundsDisc.

        If the ambient triangulation is currently oriented, then this routine
        guarantees to preserve the orientation.

        If no 2-0, 2-1 or (if requested) 3-2 moves are available, then the
        ambient triangulation will remain entirely untouched.

        Adapted from Regina's Triangulation3.simplifyToLocalMinimum() and
        SnapPea's check_for_cancellation().

        Returns:
            True if and only if the ambient triangulation was successfully
            simplified.
        """
        changed = False     # Has anything changed ever?    (Return value.)
        changedNow = True   # Did we just change something? (Loop control.)
        while changedNow:
            changedNow = False
            edgeLab = self._edgeLab()
            for edge in self._tri.edges():
                # Make sure to leave the embedded loops untouched.
                if self.isLoopEdgeIndex( edge.index() ):
                    continue

                # If requested, try a 3-2 move.
                if include32:
                    relabelling = threeTwo( edge, edgeLab )
                    if relabelling is not None:
                        changedNow = True
                        break

                # Try a 2-0 edge move.
                # This move can destroy a loop if it bounds a disc.
                relabelling = twoZero( edge, edgeLab )
                if relabelling is not None:
                    changedNow = True
                    break

                # Try a 2-1 edge move.
                # This move can destroy a loop if it bounds a disc.
                relabelling = twoOne( edge, 0, edgeLab )
                if relabelling is not None:
                    changedNow = True
                    break
                relabelling = twoOne( edge, 1, edgeLab )
                if relabelling is not None:
                    changedNow = True
                    break

            # Did we improve the triangulation? If so, then we need to update
            # the details of the embedded loops, and then check whether we can
            # make further improvements.
            if changedNow:
                changed = True
                try:
                    # If we destroyed any of the loops, then this will raise
                    # NotLoop.
                    self._setFromRelab(relabelling)
                except NotLoop:
                    # As noted above, a loop can only get destroyed if it
                    # bounds a disc.
                    raise BoundsDisc()

        # Nothing further we can do.
        return changed

    def simplifyMonotonic(self):
        """
        Uses 2-0 edge, 2-1 edge, and 3-2 moves to monotonically reduce the
        the number of tetrahedra in the ambient triangulation, while leaving
        the embedded loops untouched.

        If one of the embedded loops bounds a disc, then this routine might
        (but is not guaranteed to) raise BoundsDisc.

        If the ambient triangulation is currently oriented, then this routine
        guarantees to preserve the orientation.

        If no 2-0, 2-1 or 3-2 moves are available, then the ambient
        triangulation will remain entirely untouched.

        Adapted from Regina's Triangulation3.simplifyToLocalMinimum().

        Returns:
            True if and only if the ambient triangulation was successfully
            simplified.
        """
        # Include 3-2 moves.
        # Might raise BoundsDisc.
        return self._simplifyMonotonicImpl(True)

    def simplify(self):
        """
        Attempts to simplify the ambient triangulation, while leaving the
        embedded loops untouched.

        This routine uses minimiseVertices() and simplifyMonotonic(), in
        combination with random 4-4 moves that leave the loops untouched.

        Note that the subroutine minimiseVertices() is *not* fully implemented
        by default. Thus, subclasses that require this routine must either:
        --> override this routine; or
        --> ensure that minimiseVertices() is suitably implemented.
        In the latter case, see the documentation for minimiseVertices() for
        details on the behaviour that must be implemented.

        If one of the embedded loops bounds a disc, then this routine might
        (but is not guaranteed to) raise BoundsDisc.

        If the ambient triangulation is currently oriented, then this routine
        guarantees to preserve the orientation.

        If the attempted simplification is not successful, then the ambient
        triangulation will remain entirely untouched.

        Adapted from Regina's Triangulation3.simplify().

        Warning:
            --> Running this routine multiple times on the same triangulation
                with embedded loops may produce different results, since the
                implementation makes random decisions.

        Returns:
            True if and only if the ambient triangulation was successfully
            simplified.
        """
        RandomEngine.reseedWithHardware()

        # Work with a clone so that we can roll back changes if we are not
        # able to reduce the number of tetrahedra.
        tempTriLoops = self.clone()
        tempTri = tempTriLoops.triangulation()

        # Start by minimising vertices. This will probably increase the
        # number of tetrahedra if the number of vertices is not already
        # minimal, but hopefully the monotonic simplification saves us.
        #
        # Might raise BoundsDisc.
        tempTriLoops.minimiseVertices()
        tempTriLoops.simplifyMonotonic()

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
                if edge.index() in tempTriLoops.loopEdgeIndices():
                    # We do not want to touch any of the embedded loops.
                    continue
                for axis in range(2):
                    #NOTE Triangulation3.has44( e, ax ) was introduced in
                    #       Regina 7.4. In older versions of Regina,
                    #       equivalent functionality (checking eligibility of
                    #       the move, but not performing it) was provided by
                    #
                    #           Triangulation3.fourFourMove(
                    #               e, ax, True, False ).
                    if tempTri.has44( edge, axis ):
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
            # simplifyMonotonic() might raise BoundsDisc.
            fourFourChoice = fourFourAvailable[
                    RandomEngine.rand(availableCount) ]
            relabelling = fourFour( *fourFourChoice, tempTriLoops._edgeLab() )
            tempTriLoops._setFromRelab(relabelling)
            if tempTriLoops.simplifyMonotonic():
                # We successfully simplified!
                # Start all over again.
                fourFourAttempts = 0
                fourFourCap = 0
            else:
                fourFourAttempts += 1

        # If simplification was successful (either by reducing the number of
        # tetrahedra, or failing that by reducing the number of vertices
        # without increasing the number of tetrahedra), then sync self with
        # the now-simpler tempTriLoops.
        tri = self.triangulation()
        simplified = ( tempTri.size() < tri.size() or
                ( tempTri.size() == tri.size() and
                    tempTri.countVertices() < tri.countVertices() ) )
        if simplified:
            self.setFromLoops(tempTriLoops)
        return simplified


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

    def allowsCrush( self, surf ):
        """
        Are we able to crush this edge-ideal triangulation along the given
        normal surface?

        Letting W denote self.weight(surf), this routine returns True if and
        only if surf is of one of the following types:
        --> A 2-sphere with either W == 2 or W == 0.
        --> A disc with W == 1 and nontrivial boundary curve.
        --> A disc with W == 0.
        --> A projective plane with W == 1.

        Pre-condition:
        --> self.triangulation() is orientable.
        --> The given normal surface is embedded in self.triangulation().
        """
        surfType = SurfaceType.recognise(surf)
        weight = self.weight(surf)
        if surfType == SurfaceType.SPHERE:
            return ( weight in {0, 2} )
        elif surfType == SurfaceType.DISC:
            if weight == 1:
                return hasOnlyNonTrivialBoundaryCurves(surf)
            else:
                return ( weight == 0 )
        elif surfType == SurfaceType.RP3:
            return ( weight == 1 )
        return False

    #TODO Fix this so that it detects slope-reversing annuli.
    def splitIntoChords( self, surf ):
        """
        Returns a set containing the chords into which the given normal
        surface surf splits the union of ideal loops.

        Some or all of the ends of the returned chords will be abstractly
        joined together in pairs to indicate how all these chords would
        combine to form new loops after crushing surf. Some chords might not
        form new loops due to their ends lying in real boundary, in which
        case their ends will be left unjoined.

        Precondition:
        --> self.triangulation() is orientable.
        --> The given normal surface is embedded in self.triangulation().
        --> self.allowsCrush(surf) must be True.
        """
        chordsByLoopIndex = []
        ans = set()
        for embLoop in self:
            chords = embLoop.splitIntoChords(surf)
            chordsByLoopIndex.append(chords)
            ans.update(chords)

        # If the ideal weight is 2, then we might have two chords which are
        # joined to each other after crushing.
        #
        # From the pre-conditions, to have ideal weight 2, surf must be a
        # 2-sphere.
        incidentLoopWts = self.incidentLoopWeights(surf)
        if sum( incidentLoopWts.values() ) == 2:
            # Ideal weight 2 implies that surf is incident to exactly 2
            # ideal chords.
            incidentChords = []
            for loopInd in incidentLoopWts:
                for incChord in chordsByLoopIndex[loopInd]:
                    incidentChords.append(incChord)

            #TODO Choose one chord to translate, and get the target segments
            #       from the other chord.
            raise NotImplementedError()

        #TODO

        # Find chords (if any) that will get joined together after crushing.
        incidentLoopInds = self.incidentLoopIndices(surf)
        if len(incidentLoopInds) == 2:
            # From the preconditions, we may assume that surf is a 2-sphere.
            # We may also assume that each of the two incident loops has
            # weight 1 with respect to surf, which means that each such loop
            # is split into exactly one chord; after crushing surf, the
            # segments at the ends of the two chords will get joined
            # together, such that the chords combine to form a single new
            # loop. Note that both of the chords will be "long chords", so
            # the segments at either end of each chord are guaranteed to be
            # distinct; hence, we will have four segments that get joined to
            # each other in two pairs.
            endSegments = { loopIndex: set()
                              for loopIndex in incidentLoopInds }
            segLocations = dict()
            for loopIndex in incidentLoopInds:
                # As above, we should have exactly one chord.
                chord = chordsByLoopIndex[loopIndex].pop()
                for endNum in range(2):
                    seg = chord.endSegment(endNum)
                    endSegments[loopIndex].add(seg)
                    segLocations[seg] = ( chord, endNum )

            # Work out which two pairs of the endSegments will be joined
            # to each other after crushing surf.
            myLoopInd, yourLoopInd = endSegments.keys()
            mySeg = endSegments[myLoopInd].pop()
            yourSeg = mySeg.translateAlongSurface(
                    endSegments[yourLoopInd] )
            endSegments[yourLoopInd].remove(yourSeg)

            # Abstractly join the two pairs of endSegments together, so that
            # we can reconstruct the new embedded loops after crushing surf.
            myChord, myEndNum = segLocations[mySeg]
            yourChord, yourEndNum = segLocations[yourSeg]
            myChord.join( myEndNum, yourChord, yourEndNum )
            myChord.join( 1 - myEndNum, yourChord, 1 - yourEndNum )

        # Except for an chord which is incident to surf in the case where
        # surf is a disc, all as yet unjoined chords will join with
        # themselves to form a new loop.
        if ( SurfaceType.recognise(surf) == SurfaceType.DISC and
            len(incidentLoopInds) == 1 ):
            chordToAvoid = chordsByLoopIndex[ incidentLoopInds.pop() ].pop()
        else:
            chordToAvoid = None
        for chord in ans:
            if ( (chord != chordToAvoid) and
                (chord.joinedChord(0) is None) ):
                chord.join( 0, chord, 1 )
        return ans

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

        If no exceptions are raised, then the following are guaranteed to hold
        once this routine terminates:
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
        return super().minimiseBoundary()   # Might raise BoundsDisc.

    def _findBoundaryMove(self):
        """
        Returns details of a boundary move that simplifies the boundary of
        self.triangulation(), or None if the boundary is already minimal.

        In detail, whenever this routine returns None, the following
        conditions will all hold:
        --> Every 2-sphere boundary component will have exactly two triangles
            and three vertices.
        --> Every projective plane boundary component will have exactly two
            triangles and two vertices.
        --> Every other boundary component will have exactly one vertex.

        In the case where the boundary is not yet minimal, this routine
        guarantees to find a move that reduces the number of boundary
        triangles by two (without changing the topology of the union of ideal
        loops). The return value will be a tuple that describes this move
        using the following data:
        (0) A boundary edge e on which to perform the move.
        (1) A boolean indicating whether we need to layer across e. If this
            is True, then the move we perform will be to first layer across
            e, and then perform a close book move on the newly layered edge.
            Otherwise, the move will simply be a close book move on e.
        (2) Data for reconstructing the union of ideal loops in the new
            triangulation that results from this move, as required by the
            setFromEdgeEmbeddings() routine.

        Pre-condition:
        --> Every ideal loop is internal. Equivalently, every boundary
            triangle is disjoint from the ideal loops.
        --> Any triangular face F is incident to an ideal loop in at most one
            edge (this includes the case where multiple model edges of F are
            identified to form a single edge of some ideal loop).
            Equivalently, both of the following conditions hold:
            --- The ideal loops cannot be shortened.
            --- No triangular face forms an embedded disc that is bounded by
                one of the ideal loops.
        """
        # We can safely assume the second pre-condition because we are using
        # the default implementation for minimiseBoundary().

        # Find a boundary component that is not yet minimal.
        for bc in self._tri.boundaryComponents():
            if bc.countTriangles() <= 2 or bc.countVertices() <= 1:
                # Current boundary component is already minimal.
                continue

            # First try to find a close book move, which does not increase
            # the number of tetrahedra.
            #
            #NOTE Triangulation3.hasCloseBook(e) was introduced in
            #       Regina 7.4. In older versions of Regina, equivalent
            #       functionality (checking eligibility of the move, but
            #       not performing it) was provided by
            #
            #           Triangulation3.closeBook( e, True, False ).
            for edge in bc.edges():
                if self._tri.hasCloseBook(edge):
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

        If no exceptions are raised, then the following are guaranteed to hold
        once this routine terminates:
        --> There will be exactly one internal vertex for each ideal loop (and
            hence every ideal loop will have length 1).
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
            True if and only if this edge-ideal triangulation was changed. In
            other words, a return value of False indicates that:
            (1) the union of loops could not be shortened; and
            (2) the number of vertices in the ambient triangulation was
                already minimal to begin with.
        """
        # Can use the default implementation provided we supply an
        # implementation for _findSnapEdge().
        return super().minimiseVertices()

    def _findSnapEdge(self):
        """
        Returns details of a snap edge move that can be used to reduce the
        number of vertices in self.triangulation(), or None if the number of
        vertices is already minimal.

        In detail, whenever this routine returns None, the following
        conditions will all hold:
        --> There will be exactly one internal vertex for each ideal loop (and
            hence every ideal loop will have length 1).
        --> Every 2-sphere boundary component will have exactly two triangles
            and three vertices.
        --> Every projective plane boundary component will have exactly two
            triangles and two vertices.
        --> Every other boundary component will have exactly one vertex.

        In the case where the number of vertices is not yet minimal, this
        routine returns a tuple consisting of the following:
        (0) An edge on which a snap edge move can be performed.
        (1) Data for reconstructing the union of ideal loops in the new
            triangulation that results from this move, as required by the
            setFromEdgeEmbeddings() routine.

        Pre-condition:
        --> The ideal loops cannot be shortened.
        --> If the ambient triangulation has real boundary, then this
            boundary has already been minimised.
        """
        # We can safely assume the pre-conditions because we are using the
        # default implementation for minimiseVertices().

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

    def simplifyBasic(self):
        """
        Uses 2-0 edge and 2-1 edge moves to monotonically reduce the number
        of tetrahedra in the ambient triangulation, while leaving the ideal
        loops untouched.

        There should usually be no need to call this routine directly, since
        the functionality is subsumed by the more powerful simplify() and
        simplifyMonotonic() routines.

        If some ideal loop bounds a disc, then this routine might (but is not
        guaranteed to) raise BoundsDisc.

        If the ambient triangulation is currently oriented, then this routine
        guarantees to preserve the orientation.

        If no 2-0 or 2-1 moves are available, then the ambient triangulation
        will remain entirely untouched.

        Adapted from SnapPea's check_for_cancellation().

        Returns:
            True if and only if the ambient triangulation was successfully
            simplified.
        """
        # Do not include 3-2 moves.
        # Might raise BoundsDisc.
        return self._simplifyMonotonicImpl(False)

    def simplifyMonotonic(self):
        """
        Uses 2-0 edge, 2-1 edge, and 3-2 moves to monotonically reduce the
        number of tetrahedra in the ambient triangulation, while leaving the
        ideal loops untouched.

        There should usually be no need to call this routine directly, since
        the functionality is subsumed by the more powerful simplify() routine.

        If some ideal loop bounds a disc, then this routine might (but is not
        guaranteed to) raise BoundsDisc.

        If the ambient triangulation is currently oriented, then this routine
        guarantees to preserve the orientation.

        If no 2-0, 2-1 or 3-2 moves are available, then the ambient
        triangulation will remain entirely untouched.

        Adapted from Regina's Triangulation3.simplifyToLocalMinimum().

        Returns:
            True if and only if the ambient triangulation was successfully
            simplified.
        """
        # Just use the default implementation.
        # Might raise BoundsDisc.
        return super().simplifyMonotonic()

    def simplify(self):
        """
        Attempts to simplify the ambient triangulation, while leaving the
        ideal loops untouched.

        This routine uses minimiseVertices() and simplifyMonotonic(), in
        combination with random 4-4 moves that leave the loops untouched.

        Although this routine works very well most of the time, it can
        occasionally get stuck in a "well" that can only be escaped by
        increasing the number of tetrahedra. In such cases, it might be
        useful to try to escape using the randomise() routine.

        If one of the ideal loops bounds a disc, then this routine might (but
        is not guaranteed to) raise BoundsDisc.

        If the ambient triangulation is currently oriented, then this routine
        guarantees to preserve the orientation.

        If the attempted simplification is not successful, then the ambient
        triangulation will remain entirely untouched.

        Adapted from Regina's Triangulation3.simplify().

        Warning:
            --> Running this routine multiple times on the same edge-ideal
                triangulation may produce different results, since the
                implementation makes random decisions.

        Returns:
            True if and only if the ambient triangulation was successfully
            simplified.
        """
        # We have implemented the minimiseVertices() routine, so we can just
        # use the default implementation.
        return super().simplify()

    def randomise(self):
        """
        Attempts to randomly retriangulate this edge-ideal triangulation.

        This routine works by performing lots of random 2-3 moves, before
        attempting to simplify the ambient triangulation again (while leaving
        the ideal loops untouched). As such, this routine is often useful for
        escaping "wells" when the simplify() routine gets stuck.

        If one of the ideal loops bounds a disc, then this routine might (but
        is not guaranteed to) raise BoundsDisc.

        If the ambient triangulation is currently oriented, then this routine
        guarantees to preserve the orientation.

        Adapted from SnapPea's randomize_triangulation().
        """
        RandomEngine.reseedWithHardware()
        randomisation = 4       # Hard-coded value copied from SnapPea.
        origSize = self._tri.size()
        count = randomisation * origSize
        while count > 0:
            count -= 1

            # Attempt a random 2-3 move.
            relabelling = twoThree(
                    self._tri.triangle(
                        RandomEngine.rand( self._tri.countTriangles() ) ),
                    self._edgeLabs() )
            if relabelling is not None:
                self._setFromRelab(relabelling)

                # Try to force future random 2-3 moves to make "interesting"
                # changes.
                self.simplifyBasic()    # Might raise BoundsDisc.
                if self._tri.size() < origSize:
                    # We already succeeded in escaping the well, so we might
                    # as well terminate early.
                    break

        # Finish up by simplifying. The built-in randomness should hopefully
        # take us somewhere new.
        self.simplify()     # Might raise BoundsDisc.
        return


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
        return super().shorten()

    def minimiseBoundary(self):
        """
        Ensures that the ambient triangulation has the smallest possible
        number of boundary triangles, potentially adding tetrahedra to do
        this.

        A side-effect of calling this routine is that it will shorten the
        boundary loops if possible.

        If some boundary loop bounds a disc, then this routine might (but is
        not guaranteed to) raise BoundsDisc.

        If no exceptions are raised, then the following are guaranteed to hold
        once this routine terminates:
        --> Every boundary loop will have length 1.
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
        return super().minimiseBoundary()   # Might raise BoundsDisc.

    def _findBoundaryMove(self):
        """
        Returns details of a boundary move that simplifies the boundary of
        self.triangulation(), or None if the boundary is already minimal.

        In detail, whenever this routine returns None, the following
        conditions will all hold:
        --> Every boundary loop will have length 1.
        --> Every 2-sphere boundary component will have exactly two triangles
            and three vertices.
        --> Every projective plane boundary component will have exactly two
            triangles and two vertices.
        --> Every other boundary component will have exactly one vertex.

        In the case where the boundary is not yet minimal, this routine
        guarantees to find a move that reduces the number of boundary
        triangles by two (without changing the topology of the union of
        boundary loops). The return value will be a tuple that describes this
        move using the following data:
        (0) A boundary edge e on which to perform the move.
        (1) A boolean indicating whether we need to layer across e. If this
            is True, then the move we perform will be to first layer across
            e, and then perform a close book move on the newly layered edge.
            Otherwise, the move will simply be a close book move on e.
        (2) Data for reconstructing the union of boundary loops in the new
            triangulation that results from this move, as required by the
            setFromEdgeEmbeddings() routine.

        If some boundary loop bounds a disc, then this routine might (but is
        not guaranteed to) raise BoundsDisc.

        Pre-condition:
        --> Any boundary triangle F is incident to a boundary loop in at most
            one edge (this includes the case where multiple model edges of F
            are identified to form a single edge of some boundary loop).
            Equivalently, both of the following conditions hold:
            --- The boundary loops cannot be shortened.
            --- No boundary triangle forms an embedded disc that is bounded by
                one of the boundary loops.
        """
        # We can safely assume the pre-condition because we are using the
        # default implementation for minimiseBoundary().

        # Prioritise moves that reduce the length of the boundary loops. If
        # possible, use close book moves so that we do not introduce too many
        # new tetrahedra.
        for bloop in self:
            if len(bloop) == 1:
                continue

            # At this point, we know that bloop must lie in a boundary
            # component bc with more than one vertex (after all, bloop itself
            # has more than one vertex). Thus, bc is non-minimal except
            # possibly for the cases where bc is either a 2-sphere or a
            # projective plane. But we also have the pre-condition that bloop
            # meets each triangle of bc in at most one edge, and it is
            # straightforward to check that this rules out the minimal
            # (two-triangle) 2-spheres and projective planes.
            #
            # The upshot: bloop lies in a non-minimal boundary component.

            # We first try to find a close book move (that reduces the length
            # of bloop by 2).
            if len(bloop) > 2:
                for edge in bloop.boundaryComponent().edges():
                    # Check eligibility of close book move, but do *not*
                    # perform yet.
                    #
                    #NOTE Triangulation3.hasCloseBook(e) was introduced in
                    #       Regina 7.4. In older versions of Regina,
                    #       equivalent functionality (checking eligibility of
                    #       the move, but not performing it) was provided by
                    #
                    #           Triangulation3.closeBook( e, True, False ).
                    if not bloop.triangulation.hasCloseBook(edge):
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

        # At this point, every boundary loop has length 1. Our goal now is to
        # minimise all the boundary components without touching any of the
        # boundary loops.
        avoidEdgeIndices = self.loopEdgeIndices()
        for bc in self._tri.boundaryComponents():
            if bc.countTriangles() == 2 or bc.countVertices() == 1:
                # Already minimal.
                continue

            # First try to find a close book move.
            for edge in bc.edges():
                if edge.index() in avoidEdgeIndices:
                    continue
                #NOTE Triangulation3.hasCloseBook(e) was introduced in
                #       Regina 7.4. In older versions of Regina, equivalent
                #       functionality (checking eligibility of the move, but
                #       not performing it) was provided by
                #
                #           Triangulation3.closeBook( e, True, False ).
                if self._tri.hasCloseBook(edge):
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
                    # in the boundary loops, since we know that all the loops
                    # have length 1.
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

        # If we fell out of the boundary component loop, then all loops have
        # length 1, and all boundary components are minimal.
        return None

    def minimiseVertices(self):
        """
        Ensures that this triangulation with boundary loops has the smallest
        possible number of vertices, potentially adding tetrahedra to do this.

        A side-effect of calling this routine is that it will shorten the
        boundary loops if possible.

        If some boundary loop bounds a disc, then this routine might (but is
        not guaranteed to) raise BoundsDisc.

        If no exceptions are raised, then the following are guaranteed to hold
        once this routine terminates:
        --> Every boundary loop will have length 1.
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
        return super().minimiseVertices()

    def _findSnapEdge(self):
        """
        Returns details of a snap edge move that can be used to reduce the
        number of vertices in self.triangulation(), or None if the number of
        vertices is already minimal.

        In detail, whenever this routine returns None, the following
        conditions will all hold:
        --> Every boundary loop will have length 1.
        --> There will be no internal vertices.
        --> Every 2-sphere boundary component will have exactly two triangles
            and three vertices.
        --> Every projective plane boundary component will have exactly two
            triangles and two vertices.
        --> Every other boundary component will have exactly one vertex.

        In the case where the number of vertices is not yet minimal, this
        routine returns a tuple consisting of the following:
        (0) An edge on which a snap edge move can be performed.
        (1) Data for reconstructing the union of ideal loops in the new
            triangulation that results from this move, as required by the
            setFromEdgeEmbeddings() routine.

        Pre-condition:
        --> The boundary loops cannot be shortened.
        --> The real boundary of the ambient triangulation has already been
            minimised.
        """
        # We can safely assume the pre-conditions because we are using the
        # default implementation for minimiseVertices().

        # Find a suitable edge on which to perform a snap edge move (just
        # check whether the move is legal, do not perform yet).
        for edge in self._tri.edges():
            if snapEdge( edge, True, False ):
                # The loops should all lie entirely in the boundary, so the
                # snap edge move should never change the loops.
                return ( edge, self._edgeEmbeddingsData() )
        return

    def simplifyMonotonic(self):
        """
        Uses 2-0 edge, 2-1 edge, and 3-2 moves to monotonically reduce the
        number of tetrahedra in the ambient triangulation, while leaving the
        boundary loops untouched.

        There should usually be no need to call this routine directly, since
        the functionality is subsumed by the more powerful simplify() routine.

        If some boundary loop bounds a disc, then this routine might (but is
        not guaranteed to) raise BoundsDisc.

        If the ambient triangulation is currently oriented, then this routine
        guarantees to preserve the orientation.

        If no 2-0, 2-1 or 3-2 moves are available, then the ambient
        triangulation will remain entirely untouched.

        Adapted from Regina's Triangulation3.simplifyToLocalMinimum().

        Returns:
            True if and only if the ambient triangulation was successfully
            simplified.
        """
        # Just use the default implementation.
        # Might raise BoundsDisc.
        return super().simplifyMonotonic()

    def simplify(self):
        """
        Attempts to simplify the ambient triangulation, while leaving the
        boundary loops untouched.

        This routine uses minimiseVertices() and simplifyMonotonic(), in
        combination with random 4-4 moves (which leave the boundary loops
        untouched).

        If one of the boundary loops bounds a disc, then this routine might
        (but is not guaranteed to) raise BoundsDisc.

        If the ambient triangulation is currently oriented, then this routine
        guarantees to preserve the orientation.

        If the attempted simplification is not successful, then the ambient
        triangulation will remain entirely untouched.

        Adapted from Regina's Triangulation3.simplify().

        Warning:
            --> Running this routine multiple times on the same triangulation
                with boundary loops may produce different results, since the
                implementation makes random decisions.

        Returns:
            True if and only if the ambient triangulation was successfully
            simplified.
        """
        # We have implemented the minimiseVertices() routine, so we can just
        # use the default implementation.
        return super().simplify()
