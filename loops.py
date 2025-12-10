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
    def __init__( self, loops ):
        """
        Creates a disjoint union of the given collection of loops.

        Precondition:
        --> loops is nonempty.
        --> The elements of loops are all EmbeddedLoop objects lying inside
            the same ambient 3-manifold triangulation.
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
        cloneTri = Triangulation3(self._tri)
        cloneLoops = []
        for loop in self._loops:
            cloneLoops.append( EmbeddedLoop(
                loop._cloneImpl(cloneTri), loop.orientation() ) )
        return EmbeddedLoops(cloneLoops)

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

        The default implementation of this routine requires the helper
        routine _redirectCandidates(), which is *not* implemented by default.
        Thus, subclasses that require this routine must either:
        --> override this routine; or
        --> supply an implementation for _redirectCandidates().
        In the latter case, see the documentation for _redirectCandidates()
        for details on the behaviour that must be implemented.

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
        _shortenImpl() routine should attempt to redirect this loop.

        The EmbeddedLoops base class does not implement this routine, so
        subclasses that require this routine must provide an implementation.
        """
        raise NotImplementedError()

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

    #TODO
    pass


class IdealLoops(EmbeddedLoops):
    #TODO
    pass


class BoundaryLoops(EmbeddedLoops):
    #TODO
    pass
