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
    #TODO
    pass


class IdealLoops(EmbeddedLoops):
    #TODO
    pass


class BoundaryLoops(EmbeddedLoops):
    #TODO
    pass
