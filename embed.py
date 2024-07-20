"""
Embed a knot as an ideal loop in a triangulation of the 3-sphere.
"""
from regina import *
from loop import IdealLoop, BoundsDisc


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


def embedByFilling( knot, insertAsChild=False ):
    """
    Uses a 1/0 Dehn surgery to embed the given knot as an ideal loop in a
    triangulation of the 3-sphere.

    If insertAsChild is True and the given knot is an instance of
    PacketOfLink, then this routine will run loopPacket() on the constructed
    ideal loop, and then insert the resulting packet as a child of the given
    knot packet. This feature is switched off by default.

    This routine is mainly designed to work with *nontrivial* knots, although
    it does not check nontriviality directly. If this routine does happen to
    detect that the given knot is trivial, then it will raise BoundsDisc.

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


def reversePinch( knotComplement, packet=None ):
    """
    Builds an ideal loop representing the same knot as the given ideal
    triangulation.

    If the optional packet is supplied, then this routine will run
    loopPacket() on the constructed ideal loop, and then insert the resulting
    packet as a child of the given packet.

    This routine is mainly designed to work with *nontrivial* knots, although
    it does not check nontriviality directly. If this routine does happen to
    detect that the given knot is trivial, then it will raise BoundsDisc.

    Warning:
    --> This routine modifies the given knotComplement triangulation.
    --> This routine uses fast heuristics to attempt to construct the desired
        triangulation, and is not guaranteed to terminate.

    Returns:
        The constructed ideal loop.
    """
    # Triangulate the exterior with boundary edges appearing as the meridian
    # and longitude. The last step is not guaranteed to terminate in theory,
    # but it usually works pretty well in practice.
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
        if not loop.simplify():     # Might raise BoundsDisc.
            noSimplification += 1
    if packet is not None:
        child = loopPacket(loop)
        child.setLabel( packet.adornedLabel(
            "Embedded as edge {}".format( idealEdge.index() ) ) )
        packet.insertChildLast(child)
    return loop


def embedFromDiagram( knot, simplify=True ):
    """
    Uses a planar diagram to embed the given knot as an ideal loop in a
    triangulation of the 3-sphere.

    If simplify is True (the default), then this routine will try to simplify
    the constructed ideal loop before returning. Otherwise, if simplify is
    False, then this routine will return an ideal loop embedded in a
    triangulation of size 25 times the number of crossings in the given knot.

    This routine is mainly designed to work with *nontrivial* knots, although
    it does not check nontriviality directly. If this routine does happen to
    detect that the given knot is trivial, then it will raise BoundsDisc.

    Returns:
        The constructed ideal loop.
    """
    if knot.countComponents() > 1:
        raise ValueError( "Can only embed knots in a triangulation." )
    if knot.size() == 0:
        raise BoundsDisc()

    # Build a triangulation of the 3-sphere with a single "crossing gadget"
    # for each crossing in a planar diagram of the given knot.
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

        # Include the edge coming from the current over- or under-crossing.
        currentGadget = gadgets[ currentCrossing[0] ]
        if currentCrossing[1] % 2 == 0:
            edges.append( currentGadget.underCrossingEdge() )
        else:
            edges.append( currentGadget.overCrossingEdge() )

        # Repeat with the next crossing.
        currentCrossing = ( nextCrossing[0], (2 + nextCrossing[1]) % 4 )
        if currentCrossing == firstCrossing:
            # All done!
            loop = IdealLoop(edges)
            if simplify:
                noSimplification = 0
                while noSimplification < 2:
                    if not loop.simplify():     # Might raise BoundsDisc.
                        noSimplification += 1
            return loop
        strand = crossingToStrand[currentCrossing]


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
        # Tetrahedra 1 to 4:    Cone tetrahedra
        # Tetrahedra 5 to 8:    Outer tetrahedra
        self._tetrahedra = tri.newTetrahedra(9)

        # Join central tetrahedron to cone tetrahedra.
        for i in range(4):
            self._tetrahedra[0].join(
                    i, self._tetrahedra[1+i], Perm4(2,3) )

        # Join cone tetrahedra to each other.
        self._tetrahedra[1].join(
                3, self._tetrahedra[3], Perm4(0,3) )
        self._tetrahedra[2].join(
                2, self._tetrahedra[4], Perm4(1,2) )

        # Join outer tetrahedra to cone tetrahedra.
        self._tetrahedra[5].join(
                0, self._tetrahedra[3], Perm4(1,3,0,2) )
        self._tetrahedra[5].join(
                1, self._tetrahedra[2], Perm4(1,3,0,2) )
        self._tetrahedra[6].join(
                0, self._tetrahedra[3], Perm4(2,3,1,0) )
        self._tetrahedra[6].join(
                1, self._tetrahedra[4], Perm4(2,3,1,0) )
        self._tetrahedra[7].join(
                0, self._tetrahedra[1], Perm4(2,0,3,1) )
        self._tetrahedra[7].join(
                1, self._tetrahedra[4], Perm4(2,0,3,1) )
        self._tetrahedra[8].join(
                0, self._tetrahedra[1], Perm4(1,0,2,3) )
        self._tetrahedra[8].join(
                1, self._tetrahedra[2], Perm4(1,0,2,3) )

        # All done!
        return

    def join( self, ownStrand, other, otherStrand ):
        """
        Joins the given strand of this crossing gadget to the given strand of
        the other crossing gadget.

        Both ownStrand and otherStrand should be integers between 0 and 3,
        inclusive.
        """
        # Join own left face to other right face.
        ownLeftTet = self._tetrahedra[ 5 + ownStrand ]
        otherRightTet = other._tetrahedra[ 5 + ((otherStrand + 1) % 4) ]
        ownLeftTet.join( 3, otherRightTet, Perm4(2,3) )

        # Join own right face to other left face.
        ownRightTet = self._tetrahedra[ 5 + ((ownStrand + 1) % 4) ]
        otherLeftTet = other._tetrahedra[ 5 + otherStrand ]
        ownRightTet.join( 2, otherLeftTet, Perm4(2,3) )

        # All done!
        return

    def underCrossingEdge(self):
        """
        Returns the under-crossing edge of this crossing gadget.
        """
        return self._tetrahedra[0].edge(0,2)

    def overCrossingEdge(self):
        """
        Returns the over-crossing edge of this crossing gadget.
        """
        return self._tetrahedra[0].edge(1,3)
