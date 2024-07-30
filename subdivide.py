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
            for doomedTet in drillTet[ t.index() ].edgeSubTetrahedra(e):
                doomedIndices.add( doomedTet.index() )
    for vertexIndex in loop.vertexIndices():
        for emb in loop.triangulation().vertex(vertexIndex).embeddings():
            t = emb.tetrahedron()
            v = emb.vertex()
            for doomedTet in drillTet[ t.index() ].vertexSubTetrahedra(v):
                doomedIndices.add( doomedTet.index() )

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
    A 114-tetrahedron subdivision of a "main tetrahedron" that allows a
    neighbourhood of any subset of the main tetrahedron's 1-skeleton to be
    drilled out.
    """
    def __init__( self, tri ):
        """
        Creates a new drillable tetrahedron and inserts it into the given
        triangulation tri.
        """
        # Sub-tetrahedra 0 to 35:       Inner sub-tetrahedra
        # Sub-tetrahedra 36 to 51:      Face sub-tetrahedra
        # Sub-tetrahedra 52 to 93:      Edge sub-tetrahedra
        # Sub-tetrahedra 94 to 113:     Vertex sub-tetrahedra
        self._tet = tri.newTetrahedra(114)

        # Build groups of inner sub-tetrahedra.
        #   --> 6 groups numbered from 0 to 5, group g closest to edge g of
        #       the main tetrahedron
        #   --> each group consists of 2 central sub-tetrahedra and 4
        #       peripheral sub-tetrahedra
        #   --> self._tet[ 6*g + v ], for v in {0,1}, gives the central
        #       sub-tetrahedron closest to vertex v of edge g of the main
        #       tetrahedron
        #   --> self._tet[ 6*g + 2*v + f ], for v in {0,1} and f in {2,3},
        #       gives the peripheral sub-tetrahedron that is both closest to
        #       the triangle numbered Edge3.ordering(g)[f] and also closest
        #       to vertex v of edge g of the main tetrahedron
        for g in range(6):
            ordering = Edge3.ordering(g)

            # Join the two central sub-tetrahedra in group g.
            self._tet[ 6*g ].join( ordering[1], self._tet[ 6*g + 1 ],
                    Perm4( ordering[0], ordering[1] ) )

            # Join central sub-tetrahedra to peripheral sub-tetrahedra.
            for v in range(2):
                for f in range(2,4):
                    centralTet = self._tet[ 6*g + v ]
                    centralFace = ordering[f]
                    peripheralTet = self._tet[ 6*g + 2*v + f ]
                    gluing = Perm4( ordering[v], ordering[f] )
                    centralTet.join( centralFace, peripheralTet, gluing )

        # Build groups of face sub-tetrahedra.
        #   --> 4 groups numbered from 0 to 3, group g incident to triangle g
        #       of the main tetrahedron
        #   --> each group consists of 1 central sub-tetrahedron and 3
        #       peripheral sub-tetrahedra
        #   --> self._tet[ 36 + 5*g ] gives the central sub-tetrahedron in
        #       group g
        #   --> self._tet[ 36 + 4*g + v ], for v != g, gives the peripheral
        #       sub-tetrahedron in group g that is closest to vertex v of the
        #       main tetrahedron
        for g in range(4):
            centralTet = self._tet[ 36 + 5*g ]
            for v in range(4):
                if v == g:
                    continue
                centralFace = v
                peripheralTet = self._tet[ 36 + 4*g + v ]
                gluing = Perm4( *( {0,1,2,3} - {v,g} ) )
                centralTet.join( centralFace, peripheralTet, gluing )

        # Build groups of edge sub-tetrahedra.
        #   --> 6 groups numbered from 0 to 5, group g incident to edge g of
        #       the main tetrahedron
        #   --> each group consists of 1 central sub-tetrahedron, 2 binding
        #       sub-tetrahedra, and 4 peripheral sub-tetrahedra
        #   --> self._tet[ 52 + 7*g ] gives the central sub-tetrahedron in
        #       group g
        #   --> self._tet[ 53 + 7*g + v ], for v in {0,1}, gives the binding
        #       sub-tetrahedron that is closest to vertex v of edge g of the
        #       main tetrahedron
        #   --> self._tet[ 53 + 7*g + 2*v + f ], for v in {0,1} and f in
        #       {2,3}, gives the peripheral sub-tetrahedron that is both
        #       incident to the triangle numbered Edge3.ordering(g)[f] and
        #       also closest to vertex v of edge g of the main tetrahedron
        for g in range(6):
            ordering = Edge3.ordering(g)
            centralTet = self._tet[ 52 + 7*g ]
            for v in range(2):
                bindingTet = self._tet[ 53 + 7*g + v ]

                # Join central sub-tetrahedron to binding sub-tetrahedron.
                centralFace = ordering[1-v]
                gluing = Perm4( ordering[0], ordering[1] )
                centralTet.join( centralFace, bindingTet, gluing )

                # Join binding sub-tetrahedron to peripheral sub-tetrahedra.
                for f in range(2,4):
                    bindingFace = ordering[f]
                    peripheralTet = self._tet[ 53 + 7*g + 2*v + f ]
                    gluing = Perm4( ordering[v], ordering[f] )
                    bindingTet.join( bindingFace, peripheralTet, gluing )

        # Build groups of vertex sub-tetrahedra.
        #   --> 4 groups numbered from 0 to 3, group g incident to vertex g
        #       of the main tetrahedron
        #   --> each group consists of 1 central sub-tetrahedron, 1 apex
        #       sub-tetrahedron, and 3 base sub-tetrahedra
        #   --> self._tet[ 94 + 5*g ] gives the central sub-tetrahedron in
        #       group g
        #   --> self._tet[ 95 + 6*g ] gives the apex sub-tetrahedron in group
        #       g
        #   --> self._tet[ 95 + 5*g + v ], for v != g, gives the base
        #       sub-tetrahedron in group g that is opposite vertex v of the
        #       main tetrahedron
        flip = Perm4(2,3)
        for g in range(4):
            centralTet = self._tet[ 94 + 5*g ]
            for v in range(4):
                centralFace = flip[v]
                otherTet = self._tet[ 95 + 5*g + v ]
                centralTet.join( centralFace, otherTet, flip )

        # Join inner groups to each other.
        for e in range(5):
            myOrdering = Edge3.ordering(e)
            for v in range(2):
                for f in range(2,4):
                    ee = Edge3.faceNumber( myOrdering * Perm4(1-v,5-f) )
                    if ee < e:
                        # Only need to join from one side.
                        continue
                    yourOrdering = Edge3.ordering(ee)
                    mapping = yourOrdering.inverse() * myOrdering
                    vv = mapping[v]
                    ff = mapping[f]

                    # Join peripheral sub-tetrahedron in group e to
                    # peripheral sub-tetrahedron in group ee.
                    myTet = self._tet[ 6*e + 2*v + f ]
                    myFace = myOrdering[5-f]
                    yourTet = self._tet[ 6*ee + 2*vv + ff ]
                    gluing = Perm4( myOrdering[1-v], yourOrdering[1-vv] )
                    myTet.join( myFace, yourTet, gluing )

        # Join inner groups to face groups.
        for e in range(6):
            ordering = Edge3.ordering(e)
            for v in range(2):
                for f in range(2,4):
                    innerTet = self._tet[ 6*e + 2*v + f ]
                    innerFace = ordering[1-v]
                    faceTet = self._tet[ 36 + 4*ordering[f] + ordering[v] ]
                    gluing = Perm4( ordering[1-v], ordering[5-f] )
                    innerTet.join( innerFace, faceTet, gluing )

        # Join inner groups to edge groups.
        for g in range(6):
            ordering = Edge3.ordering(g)
            for v in range(2):
                # Join inner central sub-tetrahedra to edge binding
                # sub-tetrahedra.
                innerTet = self._tet[ 6*g + v ]
                innerFace = ordering[v]
                edgeTet = self._tet[ 53 + 7*g + v ]
                gluing = Perm4( ordering[0], ordering[1] )
                innerTet.join( innerFace, edgeTet, gluing )

                # Join inner peripheral sub-tetrahedra to edge peripheral
                # sub-tetrahedra.
                for f in range(2,4):
                    innerTet = self._tet[ 6*g + 2*v + f ]
                    innerFace = ordering[f]
                    edgeTet = self._tet[ 53 + 7*g + 2*v + f ]
                    gluing = Perm4( ordering[1-v], ordering[f] )
                    innerTet.join( innerFace, edgeTet, gluing )

        # Join edge groups to vertex groups.
        for e in range(6):
            ordering = Edge3.ordering(e)
            for v in range(2):
                for f in range(2,4):
                    edgeTet = self._tet[ 53 + 7*e + 2*v + f ]
                    edgeFace = ordering[5-f]
                    vertexTet = self._tet[ 95 + 5*ordering[v] + ordering[f] ]
                    gluing = Perm4( ordering[v], ordering[f] )
                    edgeTet.join( edgeFace, vertexTet, gluing )

        # All done!
        return

    def triangulation(self):
        """
        Returns the triangulation to which this drillable tetrahedron
        belongs.
        """
        return self._tetrahedra[0].triangulation()

    def vertexSubTetrahedra( self, v ):
        """
        Returns a list consisting of the sub-tetrahedra that form a regular
        neighbourhood of vertex v of the main tetrahedron.

        The vertex number v should be an integer between 0 and 3, inclusive.
        """
        start = 94 + 5*v
        return self._tet[start:start+5]

    def edgeSubTetrahedra( self, e ):
        """
        Returns a list consisting of the sub-tetrahedra that form a regular
        neighbourhood of the *midpoint* of edge e of the main tetrahedron.

        The edge number e should be an integer between 0 and 5, inclusive.

        To obtain a regular neighbourhood of the entire edge, it suffices to
        take the union of:
        --> the list returned by this routine; and
        --> the lists given by calling the vertexSubTetrahedra() routine on
            each of the endpoints of the edge.
        """
        start = 52 + 7*e
        return self._tet[start:start+7]

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
        # Join face sub-tetrahedra.
        for v in range(4):
            mySubTet = self._tet[ 36 + 4*myFace + v ]
            yourSubTet = other._tet[ 36 + 4*gluing[myFace] + gluing[v] ]
            mySubTet.join( myFace, yourSubTet, gluing )

        # Join edge sub-tetrahedra.
        for e in range(6):
            myOrdering = Edge3.ordering(e)
            if myFace in { myOrdering[0], myOrdering[1] }:
                continue
            f = myOrdering.inverse()[myFace]  # f in {2,3}.
            ee = Edge3.faceNumber( gluing * myOrdering )

            # Join edge central sub-tetrahedra.
            mySubTet = self._tet[ 52 + 7*e ]
            yourSubTet = other._tet[ 52 + 7*ee ]
            mySubTet.join( myFace, yourSubTet, gluing )

            # Join edge peripheral sub-tetrahedra.
            yourOrdering = Edge3.ordering(ee)
            mapping = yourOrdering.inverse() * gluing * myOrdering
            for v in range(2):
                mySubTet = self._tet[ 53 + 7*e + 2*v + f ]
                yourSubTet = other._tet[
                        53 + 7*ee + 2*mapping[v] + mapping[f] ]
                mySubTet.join( myFace, yourSubTet, gluing )

        # Join vertex sub-tetrahedra.
        for v in range(4):
            if v == myFace:
                continue

            # Join vertex apex sub-tetrahedra.
            mySubTet = self._tet[ 95 + 6*v ]
            yourSubTet = other._tet[ 95 + 6*gluing[v] ]
            mySubTet.join( myFace, yourSubTet, gluing )

            # Join vertex base sub-tetrahedra.
            mySubTet = self._tet[ 95 + 5*v + myFace ]
            yourSubTet = other._tet[ 95 + 5*gluing[v] + gluing[myFace] ]
            mySubTet.join( v, yourSubTet, gluing )

        # All done!
        return

    def linkingEdgeLocation( self, e ):
        """
        Returns a sub-tetrahedron and edge number corresponding to the link
        of edge e of the main tetrahedron.

        The edge number e should be an integer between 0 and 5, inclusive.
        """
        # The linking edge is given by edge (5-e) of an inner central
        # sub-tetrahedron.
        return ( self._tet[ 6*e ], 5-e )
