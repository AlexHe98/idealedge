"""
Auxiliary functions and classes for working with EdgeEmbedding3 objects.
"""


def edgesFromEmbeddings(edgeEmbeddings):
    """
    Converts the given container of EdgeEmbedding3 objects into a list of
    Edge3 objects.
    """
    return [ emb.tetrahedron().edge( emb.edge() ) for emb in edgeEmbeddings ]


def edgeOrientationFromEmbedding( edgeEmbedding, orientation ):
    """
    Converts the given orientation on the given edge embedding into the
    corresponding orientation on the underlying edge.

    In detail, an orientation is specified by:
    --> +1 if the edge is oriented from vertex 0 to vertex 1, or
    --> -1 if the edge is oriented from vertex 1 to vertex 0.
    A conversion is necessary because the two objects (the edge embedding and
    the underlying edge) might make different choices for which endpoints are
    labelled vertex 0 and vertex 1.
    """
    verPerm = edgeEmbedding.vertices()
    mapping = edgeEmbedding.tetrahedron().edgeMapping(
            edgeEmbedding.edge() )
    if verPerm[0] == mapping[0]:
        return orientation
    else:
        return -orientation
    return


def embeddingsFromEdgeIndices( tri, edgeIndices ):
    """
    Converts the given container of edge indices into a list of EdgeEmbedding3
    objects in the given 3-manifold triangulation.

    Note that instances of EmbeddedLoop implement all the behaviour of a
    container of edge indices, and therefore constitute valid inputs to the
    edgeIndices argument.
    """
    return [ tri.edge(ei).front() for ei in edgeIndices ]


class AbstractEdgeEmbedding:
    """
    An abstract version of Regina's EdgeEmbedding3 which references a
    tetrahedron index instead of a Tetrahedron3 object.

    Instances of this class are mutable: they support shifting the
    tetrahedron index.
    """
    def __init__( self, tetIndex, vertices ):
        """
        Creates a new abstract edge embedding containing the given data.

        Parameters
            tetIndex    The index of the tetrahedron in which the underlying
                        edge of the triangulation is contained.
            vertices    A mapping from the vertices of the underlying edge of
                        the triangulation to the corresponding vertex numbers
                        of the tetrahedron.
        """
        self._tetIndex = tetIndex
        self._vertices = vertices
        return

    @classmethod
    def makeAbstract( cls, edgeEmb ):
        """
        Creates and returns an abstract version of the given EdgeEmbedding3
        object.
        """
        return cls( edgeEmb.tetrahedron().index(), edgeEmb.vertices() )

    def tetrahedronIndex(self):
        """
        Returns the index of the tetrahedron in which the underlying edge of
        the triangulation is contained.
        """
        return self._tetIndex

    def vertices(self):
        """
        Maps vertices 0 and 1 of the underlying edge of the triangulation to
        the corresponding vertices of the tetrahedron.

        This permutation also maps 2 and 3 to the remaining vertex numbers of
        the tetrahedron.
        """
        return self._vertices

    def edge(self):
        """
        Returns the corresponding edge number of the ambient tetrahedron.

        This identifies which edge of the tetrahedron refers to the
        underlying edge of the triangulation.

        This will be between 0 and 5, inclusive.
        """
        return Edge3.faceNumber(self._vertices)

    def __eq__( self, rhs ):
        """
        Tests whether this and the given object are identical.

        Here, *identical* means that the two objects refer to the
        same-numbered tetrahedron, *and* have the same embedding permutations
        as returned by vertices().

        Since this test only examines tetrahedron and vertex numbers, it is
        meaningful to not only compare two AbstractEdgeEmbedding objects, but
        also to compare an AbstractEdgeEmbedding object on the left with an
        EdgeEmbedding3 object on the right.
        """
        if rhs is None:
            return False
        elif isinstance( rhs, AbstractEdgeEmbedding ):
            if self.tetrahedronIndex() != rhs.tetrahedronIndex():
                return False
        elif isinstance( rhs, EdgeEmbedding3 ):
            if self.tetrahedronIndex() != rhs.tetrahedron().index():
                return False
        else:
            raise TypeError( "Cannot compare equality of " +
                            "AbstractEdgeEmbedding with " +
                            "{}".format( type(rhs).__name__ ) )

        # Same tetrahedron index, so equal if and only if the vertices()
        # permutations are the same.
        return self.vertices() == rhs.vertices()

    def shiftTetrahedronIndex( self, shift ):
        """
        Shifts the tetrahedron index by the given amount.

        Under normal circumstances, we should have
            tetrahedronIndex() + shift >= 0,
        but this routine does not enforce this.
        """
        self._tetrahedronIndex += shift
        return
