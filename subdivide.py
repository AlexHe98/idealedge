"""
Tetrahedron subdivision, customised for the purpose of drilling out ideal
loops without forgetting the meridian curve.
"""
from regina import *
from loop import BoundaryLoop


def drillMeridian(loop):
    """
    Builds a BoundaryLoop representing the meridian curve obtained by
    drilling out the given IdealLoop.

    This routine might raise BoundsDisc if the meridian curve bounds a disc
    in the complement of the given loop.
    """
    # Build subdivided triangulation.
    drillTri = Triangulation3()
    drillTet = []
    for tet in loop.triangulation().tetrahedra():
        drillTet.append( DrillableTetrahedron(drillTri) )
        for face in range(4):
            adj = tet.adjacentTetrahedron(face)
            if ( adj is None or adj.index() > tet.index() or
                    ( adj.index() == tet.index() and
                        tet.adjacentFace(face) > face ) ):
                continue
            drillTet[ tet.index() ].join(
                    face, drillTet[ adj.index() ], tet.adjacentGluing(face) )

    # Find tetrahedra to drill out.
    doomedIndices = set()
    for edgeIndex in loop:
        for emb in loop.triangulation().edge(edgeIndex).embeddings():
            t = emb.tetrahedron()
            e = emb.edge()
            for doomedTet in drillTet[ t.index() ].incidentSubTetrahedra(e):
                doomedIndices.add( doomedTet.index() )
            for i in range(2):
                vertex = t.vertex( emb.vertices()[i] )
                for vemb in vertex.embeddings():
                    vti = vemb.tetrahedron().index()
                    v = vemb.vertex()
                    doomedIndices.add(
                            drillTet[vti].vertexSubTetrahedron(v).index() )

    # Find meridian loop.
    meridianLocations = []
    for emb in loop.triangulation().edge( loop[0] ).embeddings():
        tetIndex = emb.tetrahedron().index()
        edgeNum = emb.edge()
        meridianLocations.append(
                drillTet[tetIndex].linkingEdgeLocation(edgeNum) )

    # Drill.
    for i in reversed( sorted(doomedIndices) ):
        drillTri.removeTetrahedronAt(i)
    meridianEdges = [ tet.edge(edgeNum)
            for tet, edgeNum in meridianLocations ]
    meridian = BoundaryLoop(meridianEdges)
    meridian.simplify()     # Might raise BoundsDisc.
    return meridian


class DrillableTetrahedron:
    """
    A 104-tetrahedron subdivision of a tetrahedron that allows a
    neighbourhood of any subset of its 1-skeleton to be drilled out.
    """
    def __init__( self, tri ):
        """
        Creates a new drillable tetrahedron and inserts it into the given
        triangulation tri.
        """
        # Tetrahedra 0 to 19:       Inner tetrahedra
        #   --> 4 groups numbered from 0 to 3
        #   --> each group consists of 1 central tetrahedron, 1 apex
        #       tetrahedron, and 3 base tetrahedra
        #   --> self._tetrahedra[ 5*g ] gives the central tetrahedron in
        #       group g
        #   --> self._tetrahedra[ 1 + 6*g ] gives the apex tetrahedron in
        #       group g, which is the closest tetrahedron to vertex g
        #   --> self._tetrahedra[ 1 + 5*g + v ], for v != g, gives the base
        #       tetrahedron in group g that is closest to vertex v
        # Tetrahedra 20 to 35:      Face tetrahedra
        #   --> 4 groups numbered from 0 to 3, group g incident to face g
        #   --> each group consists of 1 outer tetrahedron and 3 inner
        #       tetrahedra
        #   --> self._tetrahedra[ 20 + 5*g ] gives the outer tetrahedron in
        #       group g
        #   --> self._tetrahedra[ 20 + 4*g + v ], for v != g, gives the inner
        #       tetrahedron in group g that is opposite vertex v
        # Tetrahedra 36 to 83:      Edge tetrahedra
        #   --> 6 groups numbered from 0 to 5, group g incident to edge g
        #   --> each group consists of 4 outer tetrahedra and 4 inner
        #       tetrahedra
        #   --> self._tetrahedra[ 36 + 8*g + v ] gives the outer tetrahedron
        #       closest to vertex v
        #   --> for i,j in {0,1}, self._tetrahedra[ 40 + 8*g + 2*i + j ]
        #       gives one of the two inner tetrahedra that are closer to
        #       ord[i] than ord[1-i], where ord = Edge3.ordering(g)
        # Tetrahedra 84 to 103:     Vertex tetrahedra
        #   --> 4 groups numbered from 0 to 3, group g incident to vertex g
        #   --> each group consists of 1 central tetrahedron, 1 apex
        #       tetrahedron, and 3 base tetrahedra
        #   --> self._tetrahedra[ 84 + 5*g ] gives the central tetrahedron
        #       in group g
        #   --> self._tetrahedra[ 85 + 6*g ] gives the apex tetrahedron in
        #       group g
        #   --> self._tetrahedra[ 85 + 5*g + v ], for v != g, gives the base
        #       tetrahedron in group g that is opposite vertex v
        self._tetrahedra = tri.newTetrahedra(104)
        flip = Perm4(2,3)

        # Join inner tetrahedra to each other.
        for g in range(4):
            centralTet = self._tetrahedra[ 5*g ]
            for v in range(4):
                otherTet = self._tetrahedra[ 1 + 5*g + v ]
                centralFace = flip[v]
                centralTet.join( centralFace, otherTet, flip )

        # Join face tetrahedra to each other.
         for g in range(4):
             outerTet = self._tetrahedra[ 20 + 5*g ]
             for v in range(4):
                 if v == g:
                     continue
                 innerTet = self._tetrahedra[ 20 + 4*g + v ]
                 gluing = Perm4(g,v)
                 outerTet.join( v, innerTet, gluing )

        # Join inner tetrahedra to face tetrahedra.
        for edgeNum in range(6):
            ordering = Edge3.ordering(edgeNum)
            for apex in range(2):
                g = ordering[apex]
                v = ordering[1-apex}
                innerTet = self._tetrahedra[ 1 + 5*g + v ]
                for base in range(2,4):
                    innerFace = ordering[base]
                    g = ordering[5-base]
                    v = ordering[base]
                    faceTet = self._tetrahedra[ 20 + 4*g + v ]
                    gluing =
                    innerTet.join( innerFace, faceTet, gluing )
        #TODO

        # Join edge tetrahedra to each other.
        #TODO

        # Join inner tetrahedra to edge tetrahedra.
        #TODO

        # Join face tetrahedra to edge tetrahedra.
        #TODO

        # Join vertex tetrahedra to each other.
        #TODO

        # Join edge tetrahedra to vertex tetrahedra.
        #TODO
        #TODO Reimplement.
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

#    def _vertexTet( self, v ):
#        return self._tetrahedra[ 27 + 5*v ]
#
#    def _edgeTet( self, e ):
#        return self._tetrahedra[ 5 + 3*e ]
#
#    def _faceTet( self, face, vertex ):
#        return self._tetrahedra[ 27 + 4*vertex + face ]

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
        #TODO Reimplement.
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

#    def vertexSubTetrahedron( self, vertex ):
#        """
#        Returns the sub-tetrahedron incident to the given vertex.
#
#        The given vertex number should be between 0 and 3, inclusive.
#        """
#        #TODO Reimplement.
#        return self._vertexTet(vertex)
#
#    def incidentSubTetrahedra( self, edge ):
#        """
#        Returns a list containing all of the sub-tetrahedra incident to the
#        given edge.
#
#        The given edge number should be between 0 and 5, inclusive.
#        """
#        #TODO Reimplement.
#        leclerc = []    # Nothing, just an inchident.
#
#        # Add edge tetrahedra.
#        for i in range(3):
#            leclerc.append( self._tetrahedra[ 5 + 3*edge + i ] )
#
#        # Add vertex tetrahedra.
#        ordering = Edge3.ordering(edge)
#        for i in range(2):
#            leclerc.append( self._tetrahedra[ 23 + ordering[i] ] )
#            for j in range(4):
#                if j != ordering[1-i]:
#                    leclerc.append( self._tetrahedra[
#                        27 + 4*ordering[i] + j ] )
#
#        # All done!
#        return leclerc
#
#    def linkingEdgeLocation( self, edge ):
#        """
#        Returns a sub-tetrahedron and edge number corresponding to the link
#        of the given edge in this drillable tetrahedron.
#        """
#        #TODO Reimplement.
#        edgeNum = Edge3.faceNumber(
#                Perm4(2,3) * Edge3.ordering(edge) * Perm4(2,3,0,1) )
#        return ( self._tetrahedra[0], edgeNum )
