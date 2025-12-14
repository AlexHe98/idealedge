"""
Auxiliary classes and functions for EmbeddedLoop.
"""
from regina import *


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


def edgesFromEmbeddings(edgeEmbeddings):
    """
    Converts the given container of EdgeEmbedding3 objects into a list of
    Edge3 objects.
    """
    return [ emb.tetrahedron().edge( emb.edge() ) for emb in edgeEmbeddings ]


def embeddingsFromEdgeIndices( tri, edgeIndices ):
    """
    Converts the given container of edge indices into a list of EdgeEmbedding3
    objects in the given 3-manifold triangulation.

    Note that instances of EmbeddedLoop implement all the behaviour of a
    container of edge indices, and therefore constitute valid inputs to the
    edgeIndices argument.
    """
    return [ tri.edge(ei).front() for ei in edgeIndices ]
