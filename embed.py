"""
Embed a knot as an ideal loop in a triangulation of the 3-sphere.
"""
from regina import *


class CrossingGadget:
    """
    A triangulation with edges that can be used to realise a crossing in a
    knot diagram.
    """
    # Data for finding the ideal loop.
    _underEdgeEndpoints = (0,2)
    _overEdgeEndpoints = (1,3)
    _verticalEdgeEndpoints = [ (0,3), (2,1), (3,0), (1,2) ]

    # Data for gluing gadgets together.
    _lowerLeftTet = [18,15,16,17]
    _lowerLeftVer = [
            Perm4(2,0,3,1),
            Perm4(1,2,3,0),
            Perm4(1,3,0,2),
            Perm4(2,1,0,3) ]
    _lowerRightTet = [15,16,17,18]
    _lowerRightVer = [
            Perm4(1,0,3,2),
            Perm4(1,2,0,3),
            Perm4(2,3,0,1),
            Perm4(2,1,3,0) ]
    _middleLeftTet = [5,7,9,11]
    _middleLeftVer = [
            Perm4(0,2,1,3),
            Perm4(3,1,2,0),
            Perm4(2,0,3,1),
            Perm4(1,3,0,2) ]
    _middleRightTet = [6,8,10,12]
    _middleRightVer = [
            Perm4(0,2,3,1),
            Perm4(3,1,0,2),
            Perm4(2,0,1,3),
            Perm4(1,3,2,0) ]
    _upperLeftTet = [21,22,23,24]
    _upperLeftVer = [
            Perm4(3,0,1,2),
            Perm4(1,0,2,3),
            Perm4(0,3,2,1),
            Perm4(2,3,1,0) ]
    _upperRightTet = [22,23,24,21]
    _upperRightVer = [
            Perm4(3,0,2,1),
            Perm4(1,3,2,0),
            Perm4(0,3,1,2),
            Perm4(2,0,1,3) ]

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
        self._tetrahedra[13].join( 3, self._tetrahedra[14], Perm4(1,3) )
        self._tetrahedra[13].join( 2, self._tetrahedra[15], Perm4(2,3) )
        self._tetrahedra[13].join( 0, self._tetrahedra[16], Perm4(2,3) )
        self._tetrahedra[14].join( 0, self._tetrahedra[17], Perm4(2,3) )
        self._tetrahedra[14].join( 2, self._tetrahedra[18], Perm4(2,3) )

        # Join upper tetrahedra to each other.
        self._tetrahedra[19].join( 2, self._tetrahedra[20], Perm4(0,2) )
        self._tetrahedra[19].join( 1, self._tetrahedra[21], Perm4(2,3) )
        self._tetrahedra[19].join( 3, self._tetrahedra[22], Perm4(2,3) )
        self._tetrahedra[20].join( 3, self._tetrahedra[23], Perm4(2,3) )
        self._tetrahedra[20].join( 1, self._tetrahedra[24], Perm4(2,3) )

        # All done!
        return

    def join( self, ownStrand, other, otherStrand ):
        """
        Joins the given strand of this crossing gadget to the given strand of
        the other crossing gadget.

        Both ownStrand and otherStrand should be integers between 0 and 3,
        inclusive.
        """
        # Glue lower left face.
        myLeft = self._tetrahedra[ self._lowerLeftTet[ownStrand] ]
        myLeftFace = self._lowerLeftVer[ownStrand][3]
        yourRight = other._tetrahedra[ other._lowerRightTet[otherStrand] ]
        myLeftGluing = ( other._lowerRightVer[otherStrand] *
                self._lowerLeftVer[ownStrand].inverse() )
        myLeft.join( myLeftFace, yourRight, myLeftGluing )

        # Glue lower right face.
        myRight = self._tetrahedra[ self._lowerRightTet[ownStrand] ]
        myRightFace = self._lowerRightVer[ownStrand][3]
        yourLeft = other._tetrahedra[ other._lowerLeftTet[otherStrand] ]
        myRightGluing = ( other._lowerLeftVer[otherStrand] *
                self._lowerRightVer[ownStrand].inverse() )
        myRight.join( myRightFace, yourLeft, myRightGluing )

        # Glue middle left face.
        myLeft = self._tetrahedra[ self._middleLeftTet[ownStrand] ]
        myLeftFace = self._middleLeftVer[ownStrand][3]
        yourRight = other._tetrahedra[ other._middleRightTet[otherStrand] ]
        myLeftGluing = ( other._middleRightVer[otherStrand] *
                self._middleLeftVer[ownStrand].inverse() )
        myLeft.join( myLeftFace, yourRight, myLeftGluing )

        # Glue middle right face.
        myRight = self._tetrahedra[ self._middleRightTet[ownStrand] ]
        myRightFace = self._middleRightVer[ownStrand][3]
        yourLeft = other._tetrahedra[ other._middleLeftTet[otherStrand] ]
        myRightGluing = ( other._middleLeftVer[otherStrand] *
                self._middleRightVer[ownStrand].inverse() )
        myRight.join( myRightFace, yourLeft, myRightGluing )

        # Glue upper left face.
        myLeft = self._tetrahedra[ self._upperLeftTet[ownStrand] ]
        myLeftFace = self._upperLeftVer[ownStrand][3]
        yourRight = other._tetrahedra[ other._upperRightTet[otherStrand] ]
        myLeftGluing = ( other._upperRightVer[otherStrand] *
                self._upperLeftVer[ownStrand].inverse() )
        myLeft.join( myLeftFace, yourRight, myLeftGluing )

        # Glue upper right face.
        myRight = self._tetrahedra[ self._upperRightTet[ownStrand] ]
        myRightFace = self._upperRightVer[ownStrand][3]
        yourLeft = other._tetrahedra[ other._upperLeftTet[otherStrand] ]
        myRightGluing = ( other._upperLeftVer[otherStrand] *
                self._upperRightVer[ownStrand].inverse() )
        myRight.join( myRightFace, yourLeft, myRightGluing )

        # All done!
        return

    def underCrossingEdge(self):
        """
        Returns the under-crossing edge of this crossing gadget.
        """
        return self._tetrahedra[0].edge( *self._underEdgeEndpoints )

    def overCrossingEdge(self):
        """
        Returns the over-crossing edge of this crossing gadget.
        """
        return self._tetrahedra[0].edge( *self._overEdgeEndpoints )

    def verticalEdge( self, strand ):
        """
        Returns the vertical edge at the given strand index.

        The strand index should be an integer between 0 and 3, inclusive.
        """
        return self._tetrahedra[1+strand].edge(
                *self._verticalEdgeEndpoints[strand] )
