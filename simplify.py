"""
Simplify a triangulation without destroying ideal loops.
"""
from regina import *
from moves import *


def pinchEdge(edge):
    """
    Pinches the given edge to a point.

    If the triangulation containing the given edge is oriented, then this
    operation will preserve the orientation.

    Parameters:
    --> edge    The edge that should be pinched.

    Returns:
        A dictionary mapping old edge indices to new edge indices.
    """
    #TODO
    return


def minimiseVertices( tri, idealEdgeIndices ):
    """
    ...

    Returns:
        The new set of ideal edge indices.
    """
    #TODO
    return


def monotonicSimplify( tri, idealEdgeIndices ):
    """
    Monotonically simplifies tri to some local minimum number of tetrahedra,
    while preserving the collection of ideal edges.

    This routine uses 3-2 moves, edge 2-0 moves, and 2-1 moves.

    If tri is currently oriented, then this operation will preserve the
    orientation.

    Parameters:
    --> tri                 The triangulation to simplify.
    --> idealEdgeIndices    A collection of edge indices corresponding to
                            ideal edges that should be preserved.

    Returns:
        The new set of ideal edge indices.
    """
    changed = True
    while changed:
        changed = False
        #TODO
        pass
    #TODO
    return
