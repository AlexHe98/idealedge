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
    (1) It acts as a container of EmbeddedLoop objects, which are indexed in
        an arbitrary order.
    (2) It provides methods to simplify the ambient 3-manifold triangulation,
        while ensuring that the topological embedding of the union of loops is
        always preserved.
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
        """
        cloneTri = Triangulation3(self._tri)
        cloneLoops = []
        for loop in self._loops:
            cloneLoops.append( EmbeddedLoop(
                loop._cloneImpl(cloneTri), loop.orientation() ) )
        return EmbeddedLoops(cloneLoops)

    #TODO
    pass


class IdealLoops(EmbeddedLoops):
    #TODO
    pass


class BoundaryLoops(EmbeddedLoops):
    #TODO
    pass
