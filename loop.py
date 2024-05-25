"""
Ideal loops for representing torus boundary components of a 3-manifold.
"""
from regina import *


def persistentLocation(face):
    """
    Returns a persistent identifier for the location of the given face.

    In detail, this routine returns a pair (t,f), where:
    --> t is the index of a tetrahedron meeting the given face; and
    --> f is a face number of this tetrahedron that corresponds to the given
        face.
    This identifier remains valid as long as tetrahedron indices remain
    unchanged (for example, as long as no tetrahedra are ever deleted).
    """
    emb = face.embedding(0)
    return ( emb.tetrahedron().index(), emb.face() )


class IdealLoop:
    """
    A sequence of edges representing an embedded ideal loop in a 3-manifold
    triangulation.
    """
    def __init__( self, edges ):
        """
        Creates an ideal loop from the given list of edges.

        Raises ValueError if the given list of edges does not form an
        embedded closed loop, or if the order of the edges in the given list
        does not match the order in which the edges appear in the loop.

        Pre-condition:
        --> The given list of edges is nonempty, and consists of edges that
            all belong to the same 3-manifold triangulation.
        """
        edge = edges[0]
        self._tri = edge.triangulation()
        self._edgeIndices = []
        self._vertIndices = set()
        firstVert = edge.vertex(0)
        lastVert = edge.vertex(0)
        error = False
        for edge in edges:
            self._edgeIndices.append( edge.index() )
            self._vertIndices.add( lastVert.index() )
            if edge.vertex(0) == lastVert:
                lastVert = edge.vertex(1)
            elif edge.vertex(1) == lastVert:
                lastVert = edge.vertex(0)
            else:
                error = True
                break
        if ( error or ( lastVert != firstVert ) or
                ( len( self._vertIndices() ) != len(edges) ) ):
            msg = ( "Sequence of edges does not describe an embedded " +
                    "closed loop." )
            raise ValueError(msg)
        return

    def __len__(self):
        return len( self._edgeIndices )
