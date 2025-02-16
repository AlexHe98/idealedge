"""
Ideal bouquets for representing boundary components of 3-manifold.
"""
from regina import *
from loop import EmbeddedLoop


class IdealPetal(EmbeddedLoop):
    """
    A sequence of edges representing an embedded loop that forms a petal of
    an ideal bouquet.

    A core feature of this class is that it effectively stores a list of edge
    *indices* corresponding to the edges of the ideal petal. Thus, for any
    instance P of this class, the following functionality is available:
    --> (e in P) is True if and only if P.triangulation().edge(e) is an edge
        in the petal
    --> len(P) is the number of edges in the petal
    --> for i between 0 and (len(P) - 1), inclusive, P[i] returns the index
        of the ith edge in the petal
    --> iterating through the loop yields all the edge indices in order

    An ideal petal always has a single distinguished vertex that coincides
    with the root vertex of the corresponding ideal bouquet. It is guaranteed
    that the underlying list of edge indices will order the edges so that it
    is the first and last edges that are incident to the root vertex.
    """
    #TODO
    pass
