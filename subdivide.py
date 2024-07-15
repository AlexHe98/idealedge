"""
Tetrahedron subdivision, customised for the purpose of drilling out ideal
loops without forgetting the meridian curve.
"""
from regina import *


class DrillableTetrahedron:
    """
    A 43-tetrahedron subdivision of a tetrahedron that allows a neighbourhood
    of any edge or any pair of opposite edges to be drilled out.
    """
    def __init__( self, tri ):
        """
        Creates a new drillable tetrahedron and inserts it into the given
        triangulation.
        """
        # Tetrahedra 0 to 4:    Central tetrahedra
        # Tetrahedra 5 to 22:   Edge tetrahedra
        # Tetrahedra 23 to 42:  Vertex tetrahedra
        self._tetrahedra = tri.newTetrahedra(43)
        flip = Perm4(2,3)

        # Join central tetrahedra to each other.
        for i in range(4):
            self._tetrahedra[0].join(
                    flip[i], self._tetrahedra[1+i], flip )

        # Join edge tetrahedra to each other.
        for i in range(6):
            ordering = Edge3.ordering(i)
            middleTet = self._tetrahedra[ 5 + 3*i ]
            gluing = Perm4( ordering[0], ordering[1] )
            for j in range(2):
                middleFace = ordering[j]
                endTet = self._tetrahedra[ 7 + 3*i - j ]
                middleTet.join( middleFace, endTet, gluing )

        # Join central tetrahedra to edge tetrahedra.
        for i in range(6):
            ordering = Edge3.ordering(i)
            gluing = Perm4( ordering[2], ordering[3] )
            for j in range(2):
                centralTet = self._tetrahedra[ 1 + ordering[j] ]
                centralFace = ordering[1-j]
                edgeTet = self._tetrahedra[ 6 + 3*i + j ]
                centralTet.join( centralFace, edgeTet, gluing )

        # Join vertex tetrahedra to each other.
        for i in range(4):
            innerTet = self._tetrahedra[23+i]
            for j in range(4):
                outerTet = self._tetrahedra[ 27 + 4*i + j ]
                innerTet.join(
                        flip[j], outerTet, flip )

        # Join edge tetrahedra to vertex tetrahedra.
        for i in range(6):
            ordering = Edge3.ordering(i)
            for j in range(2):
                edgeTet = self._tetrahedra[ 6 + 3*i + j ]
                for jj in range(2,4):
                    edgeFace = ordering[jj]
                    vertexTet = self._tetrahedra[
                            27 + 4*ordering[j] + ordering[jj] ]
                    gluing = Perm4( ordering[jj], ordering[5-jj] )
                    edgeTet.join( edgeFace, vertexTet, gluing )

        # All done!
        return
