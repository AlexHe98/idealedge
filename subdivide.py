"""
Tetrahedron subdivision, customised for the purpose of drilling out ideal
loops without forgetting the meridian curve.
"""
from regina import *


def subdivide(tri):
    subTri = Triangulation3()
    subTet = []
    for tet in tri.tetrahedra():
        subTet.append( DrillableTetrahedron(subTri) )
        for face in range(4):
            adj = tet.adjacentTetrahedron(face)
            if ( adj is None or adj.index() > tet.index() or
                    ( adj.index() == tet.index() and
                        tet.adjacentTriangle(face) > face ) ):
                continue
            subTet[ tet.index() ].join(
                    face, subTet[ adj.index() ], tet.adjacentGluing(face) )
    return subTri


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

    def triangulation(self):
        """
        Returns the triangulation to which this drillable tetrahedron
        belongs.
        """
        return self._tetrahedra[0].triangulation()

    def _vertexTet( self, v ):
        return self._tetrahedra[ 27 + 5*v ]

    def _edgeTet( self, e ):
        return self._tetrahedra[ 5 + 3*e ]

    def _faceTet( self, face, vertex ):
        return self._tetrahedra[ 27 + 4*vertex + face ]

    def join( self, myFace, other, gluing ):
        """
        Joins the given face of this drillable tetrahedron to some face of
        another drillable tetrahedron.

        The usage of this routine is essentially identical to the usage of
        Regina's Tetrahedron3.join() routine.

        The other drillable tetrahedron will be updated automatically (i.e.,
        you only need to call join() from one side of the gluing).

        You may join some face of this drillable tetrahedron to some
        different face of the same drillable tetrahedron (i.e., you may pass
        other == self), though you cannot join a face to itself.

        Pre-condition:
        --> This drillable tetrahedron and the other drillable tetrahedron
            both belong to the same triangulation.
        --> The given face of this drillable tetrahedron is not currently
            glued to anything.
        --> The corresponding face of the other drillable tetrahedron (i.e.,
            face number gluing[myFace] of other) is likewise not currently
            glued to anything.
        --> We are not attempting to glue a face to itself (i.e., we do not
            have both other == self and gluing[myFace] == myFace).

        Parameters:
        --> myFace  The face of this drillable tetrahedron that will be glued
                    to the other drillable tetrahedron. This face number must
                    be between 0 and 3, inclusive.
        --> other   The other drillable tetrahedron that will be glued to the
                    given face of this drillable tetrahedron.
        --> gluing  A permutation that describes how the vertices of this
                    drillable tetrahedron will map to the vertices of the
                    other drillable tetrahedron across the new gluing.
        """
        for v in range(4):
            if v == myFace:
                continue

            # Join vertex tetrahedra.
            mySubTet = self._vertexTet(v)
            yourSubTet = other._vertexTet( gluing[v] )
            mySubTet.join( myFace, yourSubTet, gluing )

            # Join face tetrahedra.
            mySubTet = self._faceTet( myFace, v )
            yourSubTet = other._faceTet( gluing[myFace], gluing[v] )
            mySubTet.join( v, yourSubTet, gluing )

        # Join edge tetrahedra.
        for e in range(6):
            ordering = Edge3.ordering(e)
            if ordering[0] == myFace or ordering[1] == myFace:
                continue
            mySubTet = self._edgeTet(e)
            yourSubTet = other._edgeTet( Edge3.faceNumber(
                gluing * ordering ) )
            mySubTet.join( myFace, yourSubTet, gluing )

        # All done!
        return
