"""
Construct triangulations of orientable Seifert fibre spaces.
"""
from regina import *


class TriSolidTorus:
    """
    An oriented triangular solid torus within a triangulation.
    """
    def __init__( self, tri ):
        """
        Constructs a triangular solid torus and inserts it into the given
        3-manifold triangulation tri.
        """
        self._tet = tri.newTetrahedra(3)
        self._tet[0].join( 1, self._tet[1], Perm4(2,3) )
        self._tet[1].join( 3, self._tet[2], Perm4(2,3) )
        self._tet[2].join( 3, self._tet[0], Perm4(1,2,3,0) )
        return

    def triangulation(self):
        """
        Returns the triangulation containing this triangular solid torus.
        """
        return self._tet[0].triangulation()
