"""
Unions of embedded loops in a 3-manifold triangulation.

The classes here provide methods for simplifying the ambient 3-manifold
triangulation, while preserving the topological embedding of the loops.
"""
from regina import *
from moves import twoThree, threeTwo, twoZero, twoOne, fourFour
from insert import snapEdge, layerOn
from loop import EmbeddedLoop, IdealLoop, BoundaryLoop
#TODO Reimplement all the simplification methods so that they can handle
#   unions of more than one embedded loop.


class EmbeddedLoops:
    """
    A disjoint union of EmbeddedLoop objects inside a single 3-manifold
    triangulation.

    This is a base class that implements common functionality for the
    IdealLoops and BoundaryLoops classes. Although this base class can be
    instantiated, the functionality it offers is much less complete than its
    aforementioned subclasses.

    This class has two core features:
    (1) It provides methods to simplify the ambient 3-manifold triangulation,
        while ensuring that the topological embedding of the union of loops is
        always preserved.
    (2) It acts as a container of EmbeddedLoop objects, which are indexed in
        an arbitrary order (but the order is kept consistent no matter how
        much the ambient triangulation is simplified). In detail, for any
        instance loops of this class:
        --> (e in loops) is True if and only if e is an EmbeddedLoop belonging
            to this union of EmbeddedLoop objects (note that equality of
            EmbeddedLoop objects is determined by their location in memory,
            and so for instance clones will not be considered equal).
        --> len(loop) is the number of EmbeddedLoop objects in this union
        --> iterating through loops yields all the EmbeddedLoop objects in
            this union, in the order in which they are indexed
        --> for i between 0 and (len(loops) - 1), inclusive, loops[i] returns
            the EmbeddedLoop at index i in this union
    """
    # Base class for loops contained in this disjoint union.
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
        Sets this to be a disjoint union of the given loops.

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
        Returns a clone of this union of embedded loops.

        The cloned union will always be embedded in a copy of
        self.triangulation()
        """
        # We use the built-in type() function to make sure that subclasses
        # will construct clones of the correct type.
        cloneTri = Triangulation3(self._tri)
        cloneLoops = []
        for loop in self._loops:
            cloneLoops.append( type(loop)(
                loop._cloneImpl(cloneTri), loop.orientation() ) )
        return type(self)(cloneLoops)

    @staticmethod
    def fromBlueprint( triEncoding, *loops ):
        """
        Constructs a collection of embedded loops using a picklable blueprint,
        as constructed by either EmbeddedLoop.blueprint() or
        EmbeddedLoops.blueprint().
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
            embLoops.append( EmbeddedLoop( edges, orientation ) )
        return EmbeddedLoops(embLoops)

    def setFromEdgeLocations( self, data ):
        """
        Sets this collection of embedded loops using the given data.

        In detail, each item in data is a pair (E, D) that specifies an
        EmbeddedLoop in a triangulation T as follows:
        --> E is a list of "edge locations". That is, each item in E is a pair
            (s, n), where s is a tetrahedron of T and n is an edge number
            (from 0 to 5, inclusive) of this tetrahedron.
        --> D specifies the orientation of the EmbeddedLoop: 
            --- It is +1 if the first edge of the loop should be oriented from
                vertex 0 to vertex 1.
            --- It is -1 if the first edge of the loop should be oriented from
                vertex 1 to vertex 0.
            --- It is 0 if this routine is allowed to choose an arbitrary
                orientation on the loop.
        All of the edges described by the edge locations must lie in the same
        triangulation.

        Raises NotLoop if any of the given lists of edge locations does not
        describe an embedded closed loop, or if the order of the edges in the
        list does not match the order in which the edges appear in the loop.

        Precondition:
        --> Each given list of edge locations is nonempty.
        --> The tetrahedra given by the first entries of each edge location
            must all belong to the same 3-manifold triangulation.
        """
        embLoops = []
        for edgeLocations, orientation in data:
            embLoops.append(
                    self._LOOP_CLASS.setFromEdgeLocations(
                        edgeLocations, orientation ) )
        self.setFromLoops(embLoops)
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
        Returns the triangulation that contains this union of embedded loops.
        """
        return self._tri

    def isBoundary(self):
        """
        Does this union of embedded loops lie entirely in the boundary of the
        ambient triangulation?
        """
        for embLoop in self:
            if not embLoop.isBoundary():
                return False
        return True

    def isInternal(self):
        """
        Does this union of embedded loops lie entirely in the interior of the
        ambient triangulation?
        """
        for embLoop in self:
            if not embLoop.isInternal():
                return False
        return True

    def countLoopEdges(self):
        """
        Counts the number of edges in this collection of embedded loops.
        """
        total = 0
        for embLoop in self:
            total += len(embLoop)
        return total

    def blueprint(self):
        """
        Returns a picklable blueprint for this collection of embedded loops.

        In detail, this routine returns a tuple
            (T, E<0>, O<0>, ..., E<L-1>, O<L-1>),
        where:
        --> L = len(self).
        --> T is Regina's tight encoding of self.triangulation().
        --> E<i> is (a copy of) the list of edge indices given by the embedded
            loop at index i in this collection.
        --> O<i> is the orientation of the embedded loop at index i of this
            collection.
        The returned blueprint can be used, via the static fromBlueprint()
        routine, to build a clone of this collection of embedded loops.
        """
        ans = [ self._tri.tightEncoding() ]
        for loop in self:
            ans.append( loop.edgeIndices() )
            ans.append( loop.orientation() )
        return tuple(ans)

    #TODO Do we need intersects() and weight()?

    def intersects( self, surf ):
        """
        Returns True if and only if this union of embedded loops has nonempty
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
        Returns the number of times this union of embedded loops intersects
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
        Shortens this collection of embedded loops by looking for triangles
        that intersect a loop in two edges, and redirecting this loop to use
        the third edge.

        If no shortening is possible, then this collection of embedded loops
        will remain entirely untouched.

        The default implementation of this routine requires the helper routine
        self._redirectCandidates(), which has the following pre-condition:
        --> For every loop in self, loop must be a subclass of EmbeddedLoop
            that implements the loop._redirectCandidates() routine.

        If one of the loops bounds a disc, then this routine might (but is not
        guaranteed to) raise BoundsDisc.

        Returns:
            True if and only if this collection of embedded loops was
            successfully shortened.
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
        _shortenImpl() routine should attempt to redirect this collection of
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
        Attempts to redirect some loop in this collection across the given
        face.

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
        collection of embedded loops entirely untouched and returns False.

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
        Ensures that the triangulation containing this collection of embedded
        loops has the smallest possible number of boundary triangles,
        potentially adding tetrahedra to do this.

        The default implementation of this routine requires the following
        helper routines:
        --> _shortenImpl()
        --> _findBoundaryMove()
        Thus, subclasses that require this routine must either:
        --> override this routine; or
        --> ensure that the aforementioned helper routines are suitably
            implemented.

        A side-effect of calling this routine is that it will shorten this
        collection of embedded loops if possible.

        If one of the loops in this collection bounds a disc, then this
        routine might (but is not guaranteed to) raise BoundsDisc.

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

        If the triangulation containing this collection of loops is currently
        oriented, then this routine guarantees to preserve the orientation.

        Adapted from Regina's Triangulation3.minimiseBoundary().

        Precondition:
        --> The ambient triangulation (i.e., self.triangulation()) is valid.

        Returns:
            True if and only if this collection of loops or the ambient
            triangulation were changed. In other words, a return value of
            False indicates that:
            (1) this collection of loops could not be shortened; and
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
            changed = True
            edge, doLayer, newEdgeIndices = moveDetails

            # Make sure we will be able to find all the edges that form the
            # loops after performing the move.
            loopLocations = []
            for indexList in newEdgeIndices:
                loopLocations.append([])
                for ei in indexList:
                    emb = self._tri.edge(ei).embedding(0)
                    loopLocations[-1].append(
                            ( emb.tetrahedron(), emb.edge() ) )

            # Perform the move, and then update this loop.
            #
            # The _findBoundaryMove() routine should already ensure legality
            # of the close book move, so no need to check before performing.
            if doLayer:
                edge = layerOn(edge).edge(5)
            self._tri.closeBook( edge, False, True )
            self.setFromEdgeLocations(loopLocations)
        return

    #TODO
    pass


class IdealLoops(EmbeddedLoops):
    """
    A disjoint union of IdealLoop objects in the interior of a single
    3-manifold triangulation.

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
        loops of this class:
        --> (e in loops) is True if and only if e is an IdealLoop belonging to
            this union of IdealLoop objects (note that equality of IdealLoop
            objects is determined by their location in memory, and so for
            instance clones will not be considered equal).
        --> len(loop) is the number of IdealLoop objects in this union
        --> iterating through loops yields all the IdealLoop objects in this
            union, in the order in which they are indexed
        --> for i between 0 and (len(loops) - 1), inclusive, loops[i] returns
            the IdealLoop at index i in this union
    """
    # Base class for loops contained in this disjoint union.
    _LOOP_CLASS = IdealLoop

    def __init__( self, loops ):
        """
        Creates a disjoint union of the given collection of ideal loops.

        Precondition:
        --> loops is nonempty.
        --> The elements of loops are all IdealLoop objects lying disjointly
            in the interior of the same ambient 3-manifold triangulation.
        """
        super().__init__(loops)
        return

    def shorten(self):
        """
        Shortens this collection of ideal loops.

        In detail, if some ideal loop meets any triangle F in exactly two
        distinct edges, then it can be shortened by replacing these two edges
        with the third edge of F.

        This routine performs such shortenings until no further shortening is
        possible. If at least one such shortening occurred, then this routine
        will return True. Otherwise, this routine will leave this collection
        of ideal loops entirely untouched, and will return False.

        If some ideal loop in this collection bounds a disc, then this routine
        might (but is not guaranteed to) raise BoundsDisc.

        Returns:
            True if and only if this collection of ideal loops was
            successfully shortened.
        """
        # IdealLoop provides an appropriate implementation of
        # _redirectCandidates(), so we can just use the default implementation
        # of shortening
        return self._shortenImpl()  # Might raise BoundsDisc.

    #TODO
    pass


class BoundaryLoops(EmbeddedLoops):
    """
    A disjoint union of BoundaryLoop objects on the boundary of a single
    3-manifold triangulation.

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
        instance loops of this class:
        --> (e in loops) is True if and only if e is a BoundaryLoop belonging
            to this union of BoundaryLoop objects (note that equality of
            BoundaryLoop objects is determined by their location in memory,
            and so for instance clones will not be considered equal).
        --> len(loop) is the number of BoundaryLoop objects in this union
        --> iterating through loops yields all the BoundaryLoop objects in
            this union, in the order in which they are indexed
        --> for i between 0 and (len(loops) - 1), inclusive, loops[i] returns
            the BoundaryLoop at index i in this union
    """
    # Base class for loops contained in this disjoint union.
    _LOOP_CLASS = BoundaryLoop

    def __init__( self, loops ):
        """
        Creates a disjoint union of the given collection of boundary loops.

        Precondition:
        --> loops is nonempty.
        --> The elements of loops are all BoundaryLoop objects lying
            disjointly on the boundary of the same ambient 3-manifold
            triangulation.
        """
        super().__init__(loops)
        return

    def shorten(self):
        """
        Shortens this collection of boundary loops.

        In detail, if some boundary loop meets any boundary triangle F in
        exactly two distinct edges, then it can be shortened by replacing
        these two edges with the third edge of F.

        This routine performs such shortenings until no further shortening is
        possible. If at least one such shortening occurred, then this routine
        will return True. Otherwise, this routine will leave this collection
        of boundary loops entirely untouched, and will return False.

        If some boundary loop in this collection bounds a disc, then this routine
        might (but is not guaranteed to) raise BoundsDisc.

        Returns:
            True if and only if this collection of ideal loops was
            successfully shortened.
        """
        # BoundaryLoop provides an appropriate implementation of
        # _redirectCandidates(), so we can just use the default implementation
        # of shortening
        return self._shortenImpl()  # Might raise BoundsDisc.

    #TODO
    pass
