"""
Oriented segments resulting from splitting edges along normal surfaces.
"""
from regina import *
from aux.quad import tetHasQuads, tetQuads


class OrientedSegment:
    """
    An oriented segment resulting from splitting an edge along a normal
    surface.

    In detail, given a normal surface S and an edge e, a segment of e is a
    connected component of e minus the intersection of S with e.

    The main purpose of this class is to facilitate tracking segments through
    the operation of crushing S. For this purpose, an important notion is that
    of a "surviving embedding" associated to a segment. To define this notion,
    consider a tetrahedron T that contains no quads of the surface S, which
    means that T contains a central cell C in the cell decomposition induced
    by cutting the triangulation along S. For each edge number n in
    {0, 1, 2, 3, 4, 5}, the cell C is incident to precisely one segment seg of
    edge n of T. We say that the embedding of edge n in T is a surviving
    embedding associated to seg. The rationale for this definition is that
    after crushing S, the segment seg will become identified with the entirety
    of the edge corresponding to this surviving embedding.
    """
    def __init__( self,
                 surface, edgeIndex, segPos, orientation, survivingEmb=None ):
        """
        Initialises an oriented segment with the given data.

        In detail, the given data should be as follows:
        --> surface         The normal surface defining the segment.
        --> edgeIndex       The index of the ambient edge (i.e., the edge
                            containing the segment).
        --> segPos          The position of the segment along the ambient
                            edge e. Specifically, the segment incident to
                            vertex 0 of e is in position 0, and then
                            subsequent segments along e are numbered in
                            increasing order. The last segment is then
                            incident to vertex 1, and is in position
                            self.edgeWeight().
        --> orientation     The orientation of the segment. This is +1 if the
                            segment is oriented from vertex 0 to vertex 1 of
                            the ambient edge, and -1 if it is oriented from
                            vertex 1 to vertex 0.
        --> survivingEmb    An optional parameter specifying a surviving
                            embedding associated to the segment.
        """
        self._surface = surface
        self._edgeIndex = edgeIndex
        self._segPos = segPos
        self._orientation = orientation
        if survivingEmb is None:
            self._checkedSurvivingEmb = False
            self._survivingEmb = None
        else:
            self._checkedSurvivingEmb = True
            self._survivingEmb = survivingEmb
        return

    def integerData(self):
        """
        Returns the triple of integer data that defines this segment.

        Specifically, the returned triple consists of the following items:
        (0) self.edgeIndex()
        (1) self.segmentPosition()
        (2) self.orientation()
        """
        return ( self._edgeIndex, self._segPos, self._orientation )

    def __eq__( self, other ):
        """
        Are self and the other segment equal?

        Returns True if and only if all of the following hold:
        --> other is also an instance of OrientedSegment,
        --> self and other both have the same surface(), and
        --> self and other both have the same integerData().
        """
        return ( isinstance( other, OrientedSegment ) and
                self.integerData() == other.integerData() and
                self.surface() == other.surface() )

    def __hash__(self):
        """
        Computes a hash for this segment.

        Note that the computation uses the associated edge index, segment
        position and orientation, but does *not* use the associated normal
        surface.
        """
        return hash( self.integerData() )

    @classmethod
    def survivors( cls, surface ):
        """
        Returns a set containing all surviving segments, with both possible
        orientations, resulting from splitting edges along the given normal
        surface.

        A surviving segment is an OrientedSegment that admits an associated
        surviving embedding.
        """
        tri = surface.triangulation()
        survivorSet = set()
        for tet in tri.tetrahedra():
            if tetHasQuads( surface, tet.index() ):
                continue

            # No quads in tet, so we will find one surviving segment along
            # each edge of tet.
            for en in range(6):
                vertexPerm = tet.edgeMapping(en)
                edgeInd = tet.edge(en).index()
                survivingSegPos = cls._survivingSegmentPosition(
                        tet, vertexPerm, surface )

                # Include both +1 and -1 orientations.
                survivorSet.add(
                        cls( surface, edgeInd, survivingSegPos,
                            1, EdgeEmbedding3(
                                tet, vertexPerm ) ) )
                survivorSet.add(
                        cls( surface, edgeInd, survivingSegPos,
                            -1, EdgeEmbedding3(
                                tet, vertexPerm * Perm4(1,0,3,2) ) ) )
        return survivorSet

    @staticmethod
    def _survivingSegmentPosition( tet, vertexPerm, surface ):
        """
        Returns the position of the surviving segment along the edge of tet
        with endpoints vertexPerm[0] and vertexPerm[1].

        Precondition:
        --> tetHasQuads( surface, tet.index() ) is False.
        """
        return surface.triangles( tet.index(), vertexPerm[0] ).pythonValue()

    def surface(self):
        """
        Returns the normal surface that defines this segment.
        """
        return self._surface

    def triangulation(self):
        """
        Returns the ambient triangulation.
        """
        return self._surface.triangulation()

    def edgeIndex(self):
        """
        Returns the index of the ambient edge.
        """
        return self._edgeIndex

    def edge(self):
        """
        Returns the ambient edge.
        """
        return self.triangulation().edge( self._edgeIndex )

    def segmentPosition(self):
        """
        Returns the position of this segment along the ambient edge.

        In detail, if e is the edge containing this segment, then the segment
        incident to vertex 0 of e is in position 0, and then subsequent
        segments along e are numbered in increasing order; the last segment is
        then incident to vertex 1, and is in position self.edgeWeight().
        """
        return self._segPos

    def orientation(self):
        """
        Returns the orientation of this segment.

        This is +1 if the segment is oriented from vertex 0 to vertex 1 of the
        ambient edge, and -1 if it is oriented from vertex 1 to vertex 0.
        """
        return self._orientation

    def survivingEmbedding(self):
        """
        Returns a surviving edge embedding associated to this segment, or None
        if no such surviving embedding exists.

        If this routine successfully returns a surviving edge embedding emb,
        then the permutation emb.vertices() will be consistent with the
        orientation of this segment, in the following sense:
        --> If self.orientation() == 1, then emb.vertices()[0] and
            emb.vertices()[1] will correspond respectively to vertices 0 and 1
            of the ambient edge.
        --> If self.orientation() == -1, then this will be flipped:
            emb.vertices()[0] and emb.vertices()[1] will correspond
            respectively to vertices 1 and 0 of the ambient edge.

        The result of this routine is cached internally, so repeated calls to
        this routine will be fast and give identical results.
        """
        if self._checkedSurvivingEmb:
            return self._survivingEmb

        # Check whether a surviving embedding exists.
        self._checkedSurvivingEmb = True
        self._survivingEmb = None   # Just in case default isn't set properly.
        for emb in self.edge().embeddings():
            tet = emb.tetrahedron()
            if tetHasQuads( self._surface, tet.index() ):
                continue
            survivingSegPos = self._survivingSegmentPosition(
                    tet, emb.vertices(), self._surface )
            if survivingSegPos != self._segPos:
                continue

            # Found a surviving embedding associated to this segment.
            if self.orientation() == 1:
                self._survivingEmb = emb
            else:   # self.orientation() == -1
                self._survivingEmb = EdgeEmbedding3(
                        tet, emb.vertices() * Perm4(1,0,3,2) )
            break
        return self._survivingEmb

    def edgeWeight(self):
        """
        Returns the weight of self.surface() on the edge containing this
        segment.

        The weight will be returned as a native Python integer.
        """
        return self._surface.edgeWeight( self._edgeIndex ).pythonValue()

    def segmentType(self):
        """
        Returns the type of this segment.

        Specifically, a segment is of type k if it has k endpoints incident to
        the normal surface that defines the segment. Thus, the type will be
        either 0, 1 or 2.
        """
        wt = self.edgeWeight()
        if wt == 0:
            return 0
        if self._segPos in { 0, wt }:
            return 1
        return 2

    def _traverseOrbit( self, targets, adjacentSegments ):
        """
        Returns a target segment that is reachable by traversing an "orbit" of
        this segment, or None if no such target is reachable.

        The "orbit" is defined and traversed using the given adjacentSegments
        function, which should take any instance seg of this class as input,
        and should yield all segments "adjacent" to seg, for whatever
        definition of "adjacent" is appropriate to define the desired orbit.

        This routine only runs membership tests on the given set of targets.
        In particular, this routine will never modify targets.
        """
        # If this segment is one of the targets, then we are already done.
        if self in targets:
            return self

        # Otherwise, we search for a target using depth-first search.
        #
        # In theory, this could be done in polynomial time using the
        # Agol-Hass-Thurston weighted orbit-counting algorithm. However, this
        # search is much easier to implement, and works very well in practice.
        tri = self.triangulation()
        segStack = [self]
        visitedSegs = set()
        while segStack:
            currentSeg = segStack.pop()
            if currentSeg in visitedSegs:
                continue

            # We haven't previously visited currentSeg, so check whether any
            # of its adjacent segments are targets, and if not check whether
            # we still need to visit these adjacent segments later.
            visitedSegs.add(currentSeg)
            for adjSeg in adjacentSegments(currentSeg):
                if adjSeg in targets:
                    return adjSeg
                elif adjSeg not in visitedSegs:
                    segStack.append(adjSeg)

        # If the search terminates without finding a target, then all of the
        # targets must be unreachable via translation along self.surface().
        return None

    def translateAlongSurface( self, targets ):
        """
        Translates this type-1 segment along self.surface(), and if possible
        returns one of the given target segments onto which this segment
        translates.

        Raises ValueError if this segment is not of type 1.

        If no target segment is reachable under such translation, then this
        routine returns None. Otherwise, this routine will return a target
        segment whose orientation is consistent with this segment's
        orientation under the translation.

        This routine only runs membership tests on the given set of targets.
        In particular, this routine will never modify targets.
        """
        if self.segmentType() != 1:
            raise ValueError(
                    "translateAlongSurface() requires type-1 segment" )
        return self._traverseOrbit( targets,
                                   _adjacentSegmentsAlongSurface )

    def translateAlongParallelCells( self, targets ):
        """
        Translates this segment along parallel cells and faces induced by
        self.surface(), and if possible returns one of the given target
        segments onto which this segment translates.

        If no target segment is reachable under such translation, then this
        routine returns None. Otherwise, this routine will return a target
        segment whose orientation is consistent with this segment's
        orientation under the translation.

        This routine only runs membership tests on the given set of targets.
        In particular, this routine will never modify targets.
        """
        return self._traverseOrbit( targets,
                                   _adjacentSegmentsAlongParallelCells )


def _adjacentSegmentsAlongSurface(seg):
    """
    Yields all oriented segments that are adjacent to seg by translation
    across elementary discs of seg.surface().

    This routine guarantees to be exhaustive, but might redundantly yield an
    adjacent segment multiple times.

    Precondition:
    --> seg.segmentType() == 1.
    """
    for emb in seg.edge().embeddings():
        teti = emb.tetrahedron().index()
        en = emb.face()
        ver = emb.vertices()

        # From the precondition, seg is a type-1 segment.
        if seg._segPos == 0:
            triangleCount = seg.surface().triangles(
                    teti, ver[0] ).pythonValue()
            if triangleCount > 0:
                for adjSeg in _adjSegsParCellsAtVert0( seg, emb ):
                    yield adjSeg

                # No more segments to yield with current emb.
                continue
        else:   # seg._segPos == seg.edgeWeight()
            triangleCount = seg.surface().triangles(
                    teti, ver[1] ).pythonValue()
            if triangleCount > 0:
                for adjSeg in _adjSegsParCellsAtVert1( seg, emb ):
                    yield adjSeg

                # No more segments to yield with current emb.
                continue

        # At this point, we know that there are no triangles at the same end
        # as the given seg.
        quadType, quadCount = tetQuads( seg.surface(), teti )
        if quadType in { en, 5 - en }:
            # This is the quad type that is disjoint from seg.edge().
            quadType, quadCount = None, 0
        if quadType is None:
            # No quads incident to seg.edge(), which means that the given
            # type-1 segment must meet a triangle at the opposite end of
            # seg.edge().
            #TODO

            # No more segments to yield with current emb.
            continue

        # 
        #TODO
        raise NotImplementedError()
    #TODO
    raise NotImplementedError()


def _adjacentSegmentsAlongParallelCells(seg):
    """
    Yields all oriented segments that are adjacent to seg by translation
    across parallel cells or faces induced by seg.surface().

    This routine guarantees to be exhaustive, but might redundantly yield an
    adjacent segment multiple times.
    """
    for emb in seg.edge().embeddings():
        teti = emb.tetrahedron().index()
        en = emb.face()
        ver = emb.vertices()

        # To locate the relevant parallel cells and faces in the current
        # tetrahedron, we need to get the normal coordinates incident to
        # seg.edge().
        f = [ seg.surface().triangles( teti, ver[i] ).pythonValue()
             for i in range(2) ]
        qType, q = tetQuads( seg.surface(), teti )
        if qType in { en, 5 - en }:
            # This is the quad type that is disjoint from seg.edge().
            qType, q = None, 0

        # Does this segment belong to a parallel cell or face in the
        # current tetrahedron? If so, then we can translate along any such
        # cells or faces to find some of the adjacent segments.
        if seg._segPos < f[0]:
            # This segment belongs to either a corner cell or a parallel
            # triangular cell at vertex ver[0].
            for adjSeg in _adjSegsParCellsAtVert0( seg, emb ):
                yield adjSeg
        elif seg._segPos > f[0] + q:
            # This segment belongs to either a corner cell or a parallel
            # triangular cell at vertex ver[1].
            for adjSeg in _adjSegsParCellsAtVert1( seg, emb ):
                yield adjSeg
        elif q > 0:
            # At this point, we have f[0] <= seg._segPos <= f[0] + q, and
            # hence this segment belongs to either a wedge cell or a
            # parallel quad cell.
            qDepth = seg._segPos - f[0]     # 0 <= qDepth <= q
            for adjSeg in _adjSegsParQuadCells(
                    seg, emb, q, qType, qDepth ):
                yield adjSeg
    return


def _adjSegsParCellsAtVert0( seg, edgeEmb ):
    """
    Yields all segments adjacent to seg via translation along a corner or
    parallel triangular cell at vertex 0 of the given edge embedding.

    Precondition:
    --> The given segment seg is incident to either a corner cell or a
        parallel triangular cell at vertex 0.
    """
    tet = edgeEmb.tetrahedron()
    ver = edgeEmb.vertices()
    for otherEnd in range(4):
        if otherEnd in { ver[0], ver[1] }:
            continue

        # The given segment seg is adjacent to a segment of the edge with
        # endpoints ver[0] and otherEnd.
        #
        #                          ver[0]
        #                             •
        #                            / \
        #               ambient edge/   \
        #                          /     \
        #                   ver[1]•       •otherEnd
        #
        enOther = Edge3.faceNumber(
                Perm4( ver[1], otherEnd ) * ver )
        eiOther = tet.edge(enOther).index()
        verOther = tet.edgeMapping(enOther)
        if verOther[0] == ver[0]:
            # Same orientation.
            yield OrientedSegment(
                    seg._surface,
                    eiOther,
                    seg._segPos,
                    seg._orientation )
        else:
            # Opposite orientation.
            wtOther = seg._surface.edgeWeight(eiOther).pythonValue()
            yield OrientedSegment(
                    seg._surface,
                    eiOther,
                    wtOther - seg._segPos,
                    -1 * seg._orientation )
    return


def _adjSegsParCellsAtVert1( seg, edgeEmb ):
    """
    Yields all segments adjacent to seg via translation along a corner or
    parallel triangular cell at vertex 1 of the given edge embedding.

    Precondition:
    --> The given segment seg is incident to either a corner cell or a
        parallel triangular cell at vertex 1.
    """
    tet = edgeEmb.tetrahedron()
    ver = edgeEmb.vertices()
    for otherEnd in range(4):
        if otherEnd in { ver[0], ver[1] }:
            continue

        # The given segment seg is adjacent to a segment of the edge with
        # endpoints ver[1] and otherEnd.
        #
        #                          ver[0]
        #                             •
        #                            /
        #               ambient edge/
        #                          /
        #                   ver[1]•-------•otherEnd
        #
        enOther = Edge3.faceNumber(
                Perm4( ver[0], otherEnd ) * ver )
        eiOther = tet.edge(enOther).index()
        verOther = tet.edgeMapping(enOther)
        if verOther[0] == ver[1]:
            # Opposite orientation.
            yield OrientedSegment(
                    seg._surface,
                    eiOther,
                    seg.edgeWeight() - seg._segPos,
                    -1 * seg._orientation )
        else:
            # Same orientation.
            wtOther = seg._surface.edgeWeight(eiOther).pythonValue()
            yield OrientedSegment(
                    seg._surface,
                    eiOther,
                    wtOther - seg.edgeWeight() + seg._segPos,
                    seg._orientation )
    return


def _adjSegsType1WedgeCells( seg, edgeEmb ):
    """
    Yields all segments adjacent to seg via translation along a quad incident
    to a wedge cell.

    Precondition:
    --> The given segment seg is of type 1, with one end incident to a quad.
    """
    #TODO
    raise NotImplementedError()


def _adjSegsParQuadCells( seg, edgeEmb, qCount, qType, qDepth ):
    """
    Yields all segments adjacent to seg via translation along a parallel cell
    or face incident to a quad.

    Precondition:
    --> The given segment seg is incident to either a wedge cell or a parallel
        quad cell.
    """
    tet = edgeEmb.tetrahedron()
    ver = edgeEmb.vertices()

    # The quads divide tet into two "sides". The edge opposite this
    # segment has endpoints lying on different sides, and we can label
    # these opposite endpoints opp[i], for i in {0,1}, so that ver[i] and
    # opp[i] lie on the same side of the quads, as shown in the diagram.
    #
    #                              ver[0]
    #                                 •
    #                                /|\
    #                   ambient edge/ | \
    #                              /__|__\
    #                             /|  |  |\
    #                      ver[1]•-|--|--|-•opp[1]
    #                             \|__|__|/
    #                              \  |  /
    #                               \ | /opposite edge
    #                                \|/
    #                                 •
    #                              opp[0]
    #
    side = [ { 0, qType + 1 } ]
    side.append( {0,1,2,3} - side[0] )
    if ver[0] not in side[0]:
        side[0], side[1] = side[1], side[0]
    side[0].remove( ver[0] )
    side[1].remove( ver[1] )
    opp = [ side[0].pop(), side[1].pop() ]

    # Find all edges of tet containing segments that are adjacent to seg.
    adjVertexPerms = []
    if qDepth < qCount:
        adjVertexPerms.append(
                Perm4( ver[0], opp[1], ver[1], opp[0] ) )
    if qDepth > 0:
        adjVertexPerms.append(
                Perm4( opp[0], ver[1], ver[0], opp[1] ) )
    if qDepth != 0 and qDepth != qCount:
        adjVertexPerms.append(
                Perm4( opp[0], opp[1], ver[0], ver[1] ) )
    for vertexPerm in adjVertexPerms:
        enAdj = Edge3.faceNumber(vertexPerm)
        eiAdj = tet.edge(enAdj).index()
        verAdj = tet.edgeMapping(enAdj)
        triangles = seg._surface.triangles(
                tet.index(), verAdj[0] ).pythonValue()
        if verAdj[0] == vertexPerm[0]:
            # Same orientation.
            yield OrientedSegment(
                    seg._surface,
                    eiAdj,
                    triangles + qDepth,
                    seg._orientation )
        else:
            # Opposite orientation.
            yield OrientedSegment(
                    seg._surface,
                    eiAdj,
                    triangles + qCount - qDepth,
                    -1 * seg._orientation )
    return
