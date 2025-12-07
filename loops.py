"""
Unions of embedded loops in a 3-manifold triangulation.

The classes here provide methods for simplifying the ambient 3-manifold
triangulation, while preserving the topological embedding of the loops.
"""
from regina import *
from moves import twoThree, threeTwo, twoZero, twoOne, fourFour
from insert import snapEdge, layerOn
from loop import EmbeddedLoop, IdealLoop, BoundaryLoop
from blueprint import triangulationBlueprint, reconstructTriangulation
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
    def fromBlueprint( size, gluings, *loops ):
        """
        Constructs a collection of embedded loops using a picklable blueprint,
        as constructed by either EmbeddedLoop.blueprint() or
        EmbeddedLoops.blueprint().
        """
        tri = reconstructTriangulation( size, gluings )
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

    def blueprint(self):
        """
        Returns a picklable blueprint for this collection of embedded loops.

        In detail, this routine returns a tuple
            (S, G, E<0>, O<0>, ..., E<L-1>, O<L-1>),
        where:
        --> L = len(self).
        --> (S,G) is the blueprint, as constructed by the
            triangulationBlueprint() routine, for tri := self.triangulation().
            Specifically:
            --- S is the size (i.e., the number of tetrahedra) of tri.
            --- G is a list such that each entry is of the form [i,f,j,p], and
                describes a gluing of tri as follows:
                --> i is a tetrahedron index of tri;
                --> f is a face number of tetrahedron i;
                --> j is the index of the tetrahedron adjacent to tetrahedron
                    i along face f; and
                --> p specifies that the gluing along face f of tetrahedron i
                    is given by the permutation Perm4.S4[p].
        --> E<i> is (a copy of) the list of edge indices given by the embedded
            loop at index i in this collection.
        --> O<i> is the orientation of the embedded loop at index i of this
            collection.
        The returned blueprint can be used, via the static fromBlueprint()
        routine, to build a clone of this collection of embedded loops.
        """
        ans = [ *triangulationBlueprint(self._tri) ]
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

    #TODO
    pass


class IdealLoops(EmbeddedLoops):
    #TODO
    pass


class BoundaryLoops(EmbeddedLoops):
    #TODO
    pass
