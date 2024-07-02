"""
Embed a knot as an ideal loop in a triangulation of the 3-sphere.
"""


class CrossingGadget:
    """
    A triangulation with edges that can be used to realise a crossing in a
    knot diagram.
    """
    def __init__( self, tri ):
        """
        Creates a new crossing gadget and inserts it into the given
        triangulation.
        """
        # Tetrahedron 0:        Central tetrahedron
        # Tetrahedra 1 to 4:    Inner middle tetrahedra
        # Tetrahedra 5 to 12:   Outer middle tetrahedra
        # Tetrahedra 13 and 14: Inner lower tetrahedra
        # Tetrahedra 15 to 18:  Outer lower tetrahedra
        # Tetrahedra 19 and 20: Inner upper tetrahedra
        # Tetrahedra 21 to 24:  Outer upper tetrahedra
        self._tetrahedra = tri.newTetrahedra(25)

        # Join central tetrahedron to inner middle tetrahedra.
        central = self._tetrahedra[0]
        for i in range(4):
            other = self._tetrahedra[1+i]
            centralFace = (i+2)%4
            central.join( centralFace, other, Perm4(2,3) )

        # Join inner middle tetrahedra to outer middle tetrahedra.
        innerFaces = [ [1,2], [3,0], [2,1], [0,3] ]
        for i in range(4):
            inner = self._tetrahedra[1+i]
            for j in range(2):
                outer = self._tetrahedra[ 5 + 2*i + j ]
                inner.join( innerFaces[i][j], outer, Perm4(2,3) )

        # Join outer middle tetrahedra to each other.
        faces = [2,3,0,1]
        for i in range(4):
            ii = (i+1)%4
            myTet = self._tetrahedra[ 6 + 2*i ]
            yourTet = self._tetrahedra[ 5 + 2*ii ]
            gluing = Perm4( faces[i], faces[ii] )
            myTet.join( faces[i], yourTet, gluing )

        # Join inner middle tetrahedra to inner lower and upper tetrahedra.
        innerFaces = [0,1,3,2]
        otherTetIndices = [19,13,20,14]
        for i in range(4):
            inner = self._tetrahedra[1+i]
            other = self._tetrahedra[ otherTetIndices[i] ]
            inner.join( innerFaces[i], other, Perm4(2,3) )

        # Join outer middle tetrahedra to outer lower and upper tetrahedra.
        otherTetIndices = [ [21,22], [15,16], [23,24], [17,18] ]
        for i in range(4):
            for j in range(2):
                middle = self._tetrahedra[ 5 + 2*i + j ]
                other = self._tetrahedra[ otherTetIndices[i][j] ]
                middle.join( i, other, Perm4(2,3) )

        # Join lower tetrahedra to each other.
        #TODO

        # Join upper tetrahedra to each other.
        #TODO
        return

    def centralTet(self):
        """
        Returns the central tetrahedron.
        """
        return self._tetrahedra[0]

    def innerMiddleTet( self, index ):
        """
        Returns the inner middle tetrahedron at the given index.

        The index should be between 0 and 3, inclusive.
        """
        return self._tetrahedra[index+1]

    def outerMiddleTet( self, index ):
        """
        Returns the outer middle tetrahedron at the given index.

        The index should be between 0 and 7, inclusive.
        """
        return self._tetrahedra[index+5]
