"""
Flexible labelling scheme for the edges of a 3-manifold triangulation.
"""
from regina import *


def edgeIndFromEmb(edgeEmbedding):
    """
    Returns the index of the underlying edge of the given EdgeEmbedding3
    object.
    """
    return edgeEmbedding.tetrahedron().edge( edgeEmbedding.edge() ).index()


class EdgeLabelling:
    """
    A labelling for (some or all of) the edges of a 3-manifold triangulation.

    Such a labelling is independent of the underlying labellings given by a
    triangulation, and is designed to help keep track of edges and their
    orientations as we modify a triangulation via elementary moves.

    In detail, an instance of this class assigns a (possibly negative) integer
    index to each edge that it tracks. Each such index is associated to an
    EdgeEmbedding3 object in the underlying triangulation.

    Since using this class creates two indepedent sets of indices for the
    edges of a triangulation, we use the following terminology to distinguish
    these two sets of edge indices:
    --> A "tracked index" is the index that is tracked by this class.
    --> An "underlying index" is the edge index given by the underlying
        labelling of the triangulation.

    For comparison with Regina's underlying edge labelling scheme:
    --> Regina only uses consecutive non-negative underlying edge indices,
        from 0 up to n - 1, where n is the total number of edges. In contrast,
        this class not only allows tracked indices to be negative, but also
        even allows the set of tracked indices to be non-consecutive.
    --> Regina assigns a unique underlying index to every edge. In contrast,
        this class may optionally leave some edges untracked, and multiple
        tracked indices may be associated to EdgeEmbedding3 objects
        corresponding to the same edge (indeed, even the EdgeEmbedding3
        objects could be equal).
    --> For each edge e, Regina can give a corresponding EdgeEmbedding3 object
        in many ways, such as by calling e.front(); these EdgeEmbedding3
        objects typically satisfy strict conditions about how they relate to
        the underlying vertex labelling in the triangulation. In contrast,
        the EdgeEmbedding3 objects tracked by this class are more flexible.

    As alluded to above, this independent and flexible labelling scheme is
    intended to help with tracking edges through elementary moves. Regina's
    implementation of moves, such as 2-3 and 3-2 moves, make no guarantees
    about how the underlying edge indices or underlying vertex labellings
    change. Therefore, tracking edges through such moves requires
    independently tracking the old edge indices and EdgeEmbedding3 objects,
    which is precisely what this class facilitates.

    Here are some specific examples of use cases, where the flexibility of
    this class is crucial:
    --> In a triangulation with n edges, performing a 2-3 move changes the
        1-skeleton by creating a new edge. This class allows any integer other
        than the ones between 0 and n - 1 (inclusive) to be assigned as the
        tracked index of the new edge. Assigning a negative tracked index
        makes it easy to distinguish the new edge from the old ones.
        Alternatively, one may also choose to simply leave the new edge
        untracked.
    --> Performing a 3-2 move changes the 1-skeleton by removing an edge. If
        we wish to preserve the edge indices for all the edges that were not
        removed, we need to allow the possibility that the set of tracked
        indices is non-consecutive.
    --> Performing a 2-0 move not only removes an edge, but also causes two
        previously distinct edges to be merged together. Here, to keep track
        of this merge, it is useful to allow two tracked indices to map to
        EdgeEmbedding3 objects corresponding to the same edge.
    """
    def __init__( self, tri, labelling=None ):
        """
        Initialises an edge labelling for the given triangulation.

        If no labelling data is supplied, then the default behaviour is to
        assign a labelling that coincides, in the following sense, with the
        underlying labelling given by tri:
        --> For every edge, the default labelling assigns the same integer
            index as the underlying labelling.
        --> The default EdgeEmbedding3 object associated to each edge e is
            just e.front().

        Otherwise, the labelling data should be a dictionary whose keys are
        integers (the tracked indices), and whose values are EdgeEmbedding3
        objects corresponding to edges in tri.

        If supplied, the given labelling data may be safely modified after
        initialising this EdgeLabelling, since this EdgeLabelling will keep
        its own private (deep) copy of the labelling data.

        However, be aware that this EdgeLabelling will keep a reference to the
        original input triangulation tri. Modifications to tri might therefore
        cause the tracked EdgeEmbedding3 objects to go out of date (if such a
        modification removes a tetrahedron that is referenced by one of the
        tracked EdgeEmbedding3 objects).
        """
        if labelling is None:
            labelling = { e.index(): e.front() for e in tri.edges() }
        self._tri = tri
        self._labelling = dict()
        for i in labelling:
            # Be extra careful not to map to None.
            if labelling[i] is not None:
                # For a deep copy, make sure to clone the embedding.
                self._labelling[i] = EdgeEmbedding3( labelling[i] )
        return

    def triangulation(self):
        """
        Returns the underlying triangulation.
        """
        return self._tri

    def cloneLabelling(self):
        """
        Returns a clone of this labelling on the same underlying
        triangulation.
        """
        # Constructor already clones the labelling data.
        return EdgeLabelling( self._tri, self._labelling )

    def clone(self):
        """
        Returns a clone of this labelling on a clone of the underlying
        triangulation.
        """
        clonedTri = Triangulation3(self._tri)
        clonedLabelling = dict()
        for i, emb in self._labelling.items():
            tet = clonedTri.tetrahedron( emb.tetrahedron().index() )
            clonedLabelling[i] = EdgeEmbedding3( tet, emb.vertices() )
        return EdgeLabelling( clonedTri, clonedLabelling )

    def __len__(self):
        """
        Returns the number of edges tracked by this labelling.
        """
        return len(self._labelling)

    def __getitem__( self, index ):
        """
        Returns the EdgeEmbedding3 object at the given index, or None if the
        given index is not tracked.
        """
        return self.get(index)

    def get( self, index ):
        """
        Returns the EdgeEmbedding3 object at the given index, or None if the
        given index is not tracked.
        """
        return self._labelling.get( index, None )

    def underlyingEdgeIndex( self, index ):
        """
        Returns the underlying index corresponding to the given tracked index,
        or None if the given index is not tracked.
        """
        emb = self[index]
        if emb is None:
            return None
        return edgeIndFromEmb(emb)

    def __setitem__( self, index, emb ):
        """
        Assigns the given EdgeEmbedding3 object emb to the given tracked
        index.
        """
        if emb is None:
            # Be extra careful not to map to None.
            self.untrack(index)
        else:
            self._labelling[index] = emb
        return

    def untrack( self, index ):
        """
        Stops tracking the edge at the given index.

        Returns True if the index status changed from tracked to untracked,
        and False if the index was already untracked to begin with.
        """
        # Status changed if and only if pop() returns an actual EdgeEmbedding3
        # object (instead of None).
        return ( self._labelling.pop( index, None ) is not None )

    def __contains__( self, index ):
        """
        Is the given index tracking an edge?
        """
        return ( index in self._labelling )

    def isTracking( self, index ):
        """
        Is the given index tracking an edge?
        """
        return ( index in self )

    def __iter__(self):
        """
        Returns an iterator over the tracked indices.
        """
        return iter(self._labelling)

    def items(self):
        """
        Returns a view of the (index, embedding) pairs of this labelling.
        """
        return self._labelling.items()

    def trackedIndices(self):
        """
        Returns a view of the indices tracked by this labelling.
        """
        return self._labelling.keys()

    def edgeEmbeddings(self):
        """
        Returns a view of the EdgeEmbedding3 objects tracked by this
        labelling.
        """
        return self._labelling.values()
