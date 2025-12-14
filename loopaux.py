"""
Auxiliary classes for EmbeddedLoop.
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


class EdgeLocation:
    """
    The location of an edge in a 3-manifold triangulation.

    In detail, such a location is specified by two pieces of data:
    --> a tetrahedron t, and
    --> an edge number e.
    The edge is given by t.edge(e). The advantage of this data is that it
    persists so long as t survives within the ambient triangulation.

    This class is a wrapper for the pair (t, e). Indexing into or iterating
    through this class is equivalent to, respectively, indexing into or
    iterating through the pair.
    """
    def __init__( self, edge, avoid=set() ):
        """
        Initialises an edge location for the given edge.

        Such a location references a particular tetrahedron of
        edge.triangulation(). If there are tetrahedra that should be avoided,
        then the indices of such tetrahedra should be provided in a container
        passed to the optional avoid argument.

        Raises ValueError if every tetrahedron incident to the given edge has
        index in the avoid container.
        """
        if not self._setEdgeLocation( edge, avoid ):
            raise ValueError(
                    "Forced to avoid every tetrahedron incident to edge." )
        return

    def update( self, avoid ):
        """
        Updates this edge location to avoid all tetrahedra with indices in the
        given avoid container.

        Returns True if and only if this update was successful.
        """
        return self._setEdgeLocation( self.edge(), avoid )

    def _setEdgeLocation( self, edge, avoid ):
        for emb in edge.embeddings():
            tet = emb.tetrahedron()
            if tet.index() not in avoid:
                self._data = ( tet, emb.edge() )
                return True
        # Falling out of the loop means that we were forced to avoid every
        # tetrahedron incident to the edge.
        return False

    def __getitem__( self, index ):
        return self._data[index]

    def __iter__(self):
        return iter( self._data )

    def triangulation(self):
        """
        Returns the triangulation in which this edge is located.
        """
        return self.tetrahedron().triangulation()

    def tetrahedron(self):
        """
        Returns the tetrahedron specifying this edge location.
        """
        return self[0]

    def edgeNumber(self):
        """
        Returns the edge number of self.tetrahedron() that specifies this edge
        location.
        """
        return self[1]

    def edge(self):
        """
        Returns the edge specified by this edge location.

        This edge is not persistent as self.triangulation() is modified. It is
        recommended to explicitly call this function each time the edge is
        required, rather than storing a reference to the edge.
        """
        return self.tetrahedron().edge( self.edgeNumber() )


#TODO Test suite.
