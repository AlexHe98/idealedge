"""
Flexible labelling scheme for the edges of a 3-manifold triangulation.
"""
from regina import *


class EdgeLabelling:
    """
    A labelling for (some or all of) the edges of a 3-manifold triangulation.

    Such a labelling is independent of the underlying labellings given by a
    triangulation, and is designed to help keep track of edges and their
    orientations as we modify a triangulation via elementary moves.

    In detail, an instance of this class assigns a (possibly negative) integer
    index to each edge that it tracks. Each such index is associated to an
    EdgeEmbedding3 object in the underlying triangulation.

    For comparison with Regina's usual edge labelling scheme:
    --> Regina only uses consecutive non-negative edge indices, from 0 up to
        n - 1, where n is the total number of edges. In contrast, this class
        not only allows negative indices, but also even allows the set of
        indices to be non-consecutive.
    --> Regina assigns a unique index to every edge. In contrast, this class
        may optionally leave some edges untracked, and multiple indices may be
        associated to EdgeEmbedding3 objects corresponding to the same edge
        (indeed, even the EdgeEmbedding3 objects could be equal).
    --> For each edge e, Regina can give a corresponding EdgeEmbedding3 object
        in many ways, such as by calling e.front(); these EdgeEmbedding3
        objects typically satisfy strict conditions about how they relate to
        the underlying vertex labelling in the triangulation. In contrast,
        the EdgeEmbedding3 objects given by this class are more flexible.

    As alluded to above, this independent and flexible labelling scheme is
    intended to help with tracking edges through elementary moves. Regina's
    implementation of moves, such as 2-3 and 3-2 moves, make no guarantees
    about how edge indices or underlying vertex labellings change. Therefore,
    tracking edges through such moves requires independently tracking the
    old edge indices and EdgeEmbedding3 objects, which is precisely what this
    class facilitates.

    Here are some specific examples of use cases, where the flexibility of
    this class is crucial:
    --> In a triangulation with n edges, performing a 2-3 move changes the
        1-skeleton by creating a new edge. This class allows any integer other
        than the ones between 0 and n - 1 (inclusive) to be assigned as the
        index of the new edge. Assigning a negative index makes it easy to
        distinguish the new edge from the old ones. Alternatively, one may
        also choose to simply leave the new edge untracked.
    --> Performing a 3-2 move changes the 1-skeleton by removing an edge. If
        we wish to preserve the edge indices for all the edges that were not
        removed, we need to allow the possibility that the set of edge indices
        is non-consecutive.
    --> Performing a 2-0 move not only removes an edge, but also causes two
        previously distinct edges to be merged together. Here, to keep track
        of this merge, it is useful to allow two edge indices to map to
        EdgeEmbedding3 objects corresponding to the same edge.
    """
    def __init__( self, tri ):
        """
        Initialises the default edge labelling for the given triangulation.

        The default coincides with the underlying labelling as follows:
        --> For every edge, the default labelling assigns the same integer
            index as the underlying labelling.
        --> The default EdgeEmbedding3 object associated to each edge e is
            just e.front().
        """
        self._tri = tri
        self._labelling = { e.index(): e.front() for e in tri.edges() }
        return

    #TODO
