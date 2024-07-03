"""
Embed a knot as an ideal loop in a triangulation of the 3-sphere.
"""
from regina import *
from loop import IdealLoop


def loopPacket(loop):
    """
    Returns a packet of the triangulation containing the given loop, with an
    ideal triangulation of the drilled 3-manifold as a child.
    """
    drilled = PacketOfTriangulation3( loop.drill() )
    drilled.setLabel( "Drilled: {}".format( drilled.isoSig() ) )
    packet = PacketOfTriangulation3( loop.triangulation() )
    packet.insertChildLast(drilled)
    return packet


def reversePinch( knotComplement, packet=None ):
    """
    Builds an ideal loop representing the same knot as the given ideal
    triangulation.

    Warning:
    --> This routine modifies the given knotComplement triangulation.
    --> This routine uses fast heuristics to attempt to construct the desired
        triangulation, and is not guaranteed to terminate.

    Returns:
        The constructed ideal loop.
    """
    # Triangulate the exterior with boundary edges appearing as the meridian
    # and longitude. The last step is not guaranteed to terminate in theory,
    # but it should be fine in practice.
    knotComplement.minimiseVertices()
    knotComplement.intelligentSimplify()
    knotComplement.intelligentSimplify()
    knotComplement.idealToFinite()
    noSimplification = 0
    while noSimplification < 2:
        if not knotComplement.intelligentSimplify():
            noSimplification += 1
    mer, lon = knotComplement.meridianLongitude()

    # Get a tetrahedron index and edge number for the longitude, so that we
    # can remember its location after closing up the boundary.
    emb = lon.embedding(0)
    tet = emb.tetrahedron()
    edgeNum = emb.face()

    # Close up the boundary and build the IdealLoop.
    layer = knotComplement.layerOn(mer)
    layer.join( 0, layer, Perm4(0,1) )
    idealEdge = tet.edge(edgeNum)
    loop = IdealLoop( [idealEdge] )
    noSimplification = 0
    while noSimplification < 2:
        if not loop.simplify():
            noSimplification += 1
    if packet is not None:
        child = loopPacket(loop)
        child.setLabel( packet.adornedLabel(
            "Embedded as edge {}".format( idealEdge.index() ) ) )
        packet.insertChildLast(child)
    return loop


def embedByFilling( knot, insertAsChild=False ):
    """
    Uses a 1/0 Dehn surgery to embed the given knot as an ideal loop in a
    triangulation of the 3-sphere.

    Warning:
    --> This routine uses fast heuristics to attempt to construct the desired
        triangulation, and is not guaranteed to terminate. Use the
        embedFromDiagram() routine to build the ideal loop using an algorithm
        that guarantees to terminate.

    Returns:
        The constructed ideal loop.
    """
    if knot.countComponents() > 1:
        raise ValueError( "Can only embed knots in a triangulation." )
    if insertAsChild and isinstance( knot, PacketOfLink ):
        packet = knot
    else:
        packet = None
    knotComplement = knot.complement()
    return reversePinch( knotComplement, packet )


def embedFromDiagram( knot, simplify=True ):
    """
    Uses a planar diagram to embed the given knot as an ideal loop in a
    triangulation of the 3-sphere.

    If simplify is True (the default), then this routine will try to simplify
    the constructed ideal loop before returning. Otherwise, if simplify is
    False, then this routine will return an ideal loop embedded in a
    triangulation of size 25 times the number of crossings in the given knot.

    In the special case where the given knot has no crossings, this routine
    simply returns an unknotted ideal loop embedded in a triangulation with
    just one tetrahedron.

    Returns:
        The constructed ideal loop.
    """
    if knot.countComponents() > 1:
        raise ValueError( "Can only embed knots in a triangulation." )

    # If no crossings, use the unknotted edge in the 1-tetrahedron 3-sphere.
    if knot.size() == 0:
        tri = Triangulation3.fromIsoSig( "bkaagj" )
        return IdealLoop( [ tri.edge(0) ] )

    # Build the triangulation.
    tri = Triangulation3()
    crossings = knot.pdData()
    gadgets = []
    strandToCrossing = dict()
    crossingToStrand = dict()
    for i, crossing in enumerate(crossings):
        gadgets.append( CrossingGadget(tri) )
        for j, strand in enumerate(crossing):
            if strand in strandToCrossing:
                # We have already seen the current strand from the other end,
                # which means that we need to join the two corresponding
                # crossing gadgets.
                ii, jj = strandToCrossing[strand][0]
                gadgets[i].join( j, gadgets[ii], jj )
            else:
                strandToCrossing[strand] = []
            strandToCrossing[strand].append( (i,j) )
            crossingToStrand[ (i,j) ] = strand

    # Find the edges that form the ideal loop.
    edges = []
    firstCrossing = strandToCrossing[1][0]
    currentCrossing = strandToCrossing[1][0]
    strand = 1
    while True:
        # Find the next crossing.
        if currentCrossing == strandToCrossing[strand][0]:
            nextCrossing = strandToCrossing[strand][1]
        else:
            nextCrossing = strandToCrossing[strand][0]

        # Include the edge coming from the current crossing, and check
        # whether an extra edge is needed to ensure that the current edge is
        # connected to the edge that will come from the next crossing.
        currentGadget = gadgets[ currentCrossing[0] ]
        if currentCrossing[1] % 2 == 0:
            edges.append( currentGadget.underCrossingEdge() )
        else:
            edges.append( currentGadget.overCrossingEdge() )
        if currentCrossing[1] % 2 != nextCrossing[1] % 2:
            # Either current edge is over-crossing and next edge is
            # under-crossing, or vice versa. We use an extra vertical edge
            # to connect the current edge to the next edge.
            edges.append( currentGadget.verticalEdge( currentCrossing[1] ) )

        # Repeat with the next crossing.
        currentCrossing = ( nextCrossing[0], (2 + nextCrossing[1]) % 4 )
        if currentCrossing == firstCrossing:
            # All done!
            loop = IdealLoop(edges)
            if simplify:
                noSimplification = 0
                while noSimplification < 2:
                    if not loop.simplify():
                        noSimplification += 1
            return loop
        strand = crossingToStrand[currentCrossing]


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
