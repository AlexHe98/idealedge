"""
Oriented segments resulting from splitting edges along normal surfaces.
"""
from enum import Enum, auto
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

        See the survivingEmbedding() routine for the conditions that must be
        satisfied by the optional survivingEmb argument.
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

    def _traverseOrbit( self, targets, adjacentSegments, *args ):
        """
        Returns a target segment that is reachable by traversing an "orbit" of
        this segment, or None if no such target is reachable.

        The "orbit" is defined and traversed using the given
        adjacentSegments() function. The adjacentSegments() function should
        take as input any instance seg of this class, followed by all
        additional args supplied to this routine. Note that for a single call
        to this routine, adjacentSegments() may be called with many instances
        of this class as the first parameter, but then all subsequent
        parameters will be fixed to whatever the supplied args happen to be.
        With such input, the adjacentSegments() function should yield all
        segments "adjacent" to seg, for whatever definition of "adjacent" is
        appropriate to transitively define the desired orbit.

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
            for adjSeg in adjacentSegments( currentSeg, *args ):
                if adjSeg in targets:
                    return adjSeg
                elif adjSeg not in visitedSegs:
                    segStack.append(adjSeg)

        # If the search terminates without finding a target, then all of the
        # targets must be unreachable via translation along self.surface().
        return None

    def _getEndIncidences( self, edgeEmb ):
        r"""
        Classifies the objects incident to the endpoints of this segment,
        relative to the given edge embedding.

        In detail, the edge embedding should describe how self.edge() is
        embedded in some tetrahedron. The return value will then be a
        2-element list R structured as follows. For each i in {0, 1}, consider
        the end of seg that is closer to vertex i of the given edge embedding.
        R[i] will be a pair P that specifies the object incident to this end
        as follows:
        --> If the end is incident to a vertex v, then P will be
                ( _SegmentEndIncidence.VERTEX, j ),
            where j in {0, 1} is the vertex number assigned to v by edgeEmb.
        --> If the end is incident to a normal triangle t, then P will be
                ( _SegmentEndIncidence.TRIANGLE, j ),
            where j in {0, 1} is vertex number assigned by edgeEmb to the
            vertex that is cut off by t.
        --> If the end is incident to a normal quad q, then P will be
                ( _SegmentEndIncident.QUAD, k ),
            where k in {0, 1, 2} indicates that q separates the edges number k
            and 5 - k.
        """
        teti = edgeEmb.tetrahedron().index()
        en = edgeEmb.edge()
        ver = edgeEmb.vertices()
        ans = []

        # The answer depends on where this segment sits relative to the
        # elementary discs.
        triangleCount = [
                self.surface().triangles( teti, ver[i] ).pythonValue()
                for i in range(2) ]
        quadType, quadCount = tetQuads( self.surface(), teti )
        if ( quadType is None ) or ( quadType in { en, 5 - en } ):
            # Even if there is a quad, it won't incident to self.edge(), so we
            # don't want to count it.
            quadCount = 0

        # First incidence.
        if self._segPos == 0:
            ans.append( ( _SegmentEndIncidence.VERTEX, 0 ) )
        elif self._segPos <= triangleCount[0]:
            ans.append( ( _SegmentEndIncidence.TRIANGLE, 0 ) )
        elif self._segPos <= triangleCount[0] + quadCount:
            ans.append( ( _SegmentEndIncidence.QUAD, quadType ) )
        else:
            ans.append( ( _SegmentEndIncidence.TRIANGLE, 1 ) )

        # Second incidence.
        if self._segPos < triangleCount[0]:
            ans.append( ( _SegmentEndIncidence.TRIANGLE, 0 ) )
        elif self._segPos < triangleCount[0] + quadCount:
            ans.append( ( _SegmentEndIncidence.QUAD, quadType ) )
        elif self._segPos < self.edgeWeight():
            ans.append( ( _SegmentEndIncidence.TRIANGLE, 1 ) )
        else:
            ans.append( ( _SegmentEndIncidence.VERTEX, 1 ) )

        # All done!
        return ans

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

        # Mark the end of this segment that is incident to self.surface().
        if self._segPos == 0:
            if self._orientation == 1:
                markedEnd = 1
            else:   # self._orientation == -1
                markedEnd = 0
        else:   # self._segPos == self.edgeWeight()
            if self._orientation == 1:
                markedEnd = 0
            else:   # self._orientation == -1
                markedEnd = 1
        return self._traverseOrbit( targets,
                                   _adjacentSegmentsAlongSurface,
                                   markedEnd )

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


class _SegmentEndIncidence(Enum):
    """
    Classification of the object incident to an endpoint of a segment.

    This is either:
    --> a vertex of the ambient edge,
    --> a normal triangle touching the ambient edge, or
    --> a normal quad touching the ambient edge.
    """
    VERTEX = auto()
    TRIANGLE = auto()
    QUAD = auto()
    pass


def _adjacentSegmentsAlongSurface( seg, markedEnd ):
    """
    Yields all oriented segments that are adjacent to seg by translating the
    markedEnd across elementary discs of seg.surface().

    The markedEnd specifies one of the two endpoints of seg relative to its
    orientation. Specifically:
    --> if seg is oriented away from the marked endpoint, then markedEnd
        should be 0; and
    --> if seg is oriented towards the marked endpoint, then markedEnd should
        be 1.

    This routine guarantees to be exhaustive, but might redundantly yield an
    adjacent segment multiple times.

    This routine raises ValueError if the markedEnd of seg is not incident to
    seg.surface().
    """
    for edgeEmb in seg.edge().embeddings():
        ver = edgeEmb.vertices()
        endIncidences = seg._getEndIncidences(edgeEmb)

        # The endIncidences pair labels the ends 0 and 1 according to the
        # labelling of seg.edge() given by the ver permutation. In contrast,
        # the given markedEnd labels the ends 0 and 1 according to the
        # orientation of the segment, which might differ from the ver
        # labelling.
        if seg._orientation == 1:
            markedIncidence = endIncidences[markedEnd]
        else:   # seg._orientation == -1
            markedIncidence = endIncidences[ 1 - markedEnd ]

        # The adjacent segments depend on whether the markedEnd is incident to
        # a triangle or a quad.
        if markedIncidence[0] == _SegmentEndIncidence.TRIANGLE:
            if markedIncidence[1] == 0:
                for adjSeg in _adjSegsAlongTriangleAtVert0( seg, edgeEmb ):
                    yield adjSeg
            else:   # markedIncidence[1] == 1
                for adjSeg in _adjSegsAlongTriangleAtVert1( seg, edgeEmb ):
                    yield adjSeg
        elif markedIncidence[0] == _SegmentEndIncidence.QUAD:
            includeNonParallel = True
            for adjSeg in _adjSegsAlongQuad(
                    seg, edgeEmb, markedIncident[1], includeNonParallel ):
                yield adjSeg
        else:
            raise ValueError(
                    "The markedEnd should be incident to seg.surface()" )
    return


def _adjacentSegmentsAlongParallelCells(seg):
    """
    Yields all oriented segments that are adjacent to seg by translation
    across parallel cells or faces induced by seg.surface().

    This routine guarantees to be exhaustive, but might redundantly yield an
    adjacent segment multiple times.
    """
    for edgeEmb in seg.edge().embeddings():
        endIncidences = seg._getEndIncidences(edgeEmb)

        # Does this segment belong to a parallel cell or face in the
        # current tetrahedron? If so, then we can translate along any such
        # cells or faces to find some of the adjacent segments.
        includeNonParallel = False
        if endIncidences[1] == ( _SegmentEndIncidence.TRIANGLE, 0 ):
            # This segment belongs to either a corner cell or a parallel
            # triangular cell at vertex edgeEmb.vertices()[0].
            for adjSeg in _adjSegsAlongTriangleAtVert0(
                    seg, edgeEmb ):
                yield adjSeg
        elif endIncidences[0] == ( _SegmentEndIncidence.TRIANGLE, 1 ):
            # This segment belongs to either a corner cell or a parallel
            # triangular cell at vertex edgeEmb.vertices()[1].
            for adjSeg in _adjSegsAlongTriangleAtVert1(
                    seg, edgeEmb ):
                yield adjSeg
        elif endIncidences[0][0] == _SegmentEndIncidence.QUAD:
            # This segment is incident to a quad.
            for adjSeg in _adjSegsAlongQuad(
                    seg, edgeEmb, endIncidences[0][1], includeNonParallel ):
                yield adjSeg
        elif endIncidences[1][0] == _SegmentEndIncidence.QUAD:
            # Again, this segment is incident to a quad.
            for adjSeg in _adjSegsAlongQuad(
                    seg, edgeEmb, endIncidences[1][1], includeNonParallel ):
                yield adjSeg
    return


def _adjSegsAlongTriangleAtVert0( seg, edgeEmb ):
    """
    Yields all segments adjacent to seg via translation along a normal
    triangle at vertex 0 of the given edge embedding.

    Precondition:
    --> The given segment seg is incident to a normal triangle at vertex 0.
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


def _adjSegsAlongTriangleAtVert1( seg, edgeEmb ):
    """
    Yields all segments adjacent to seg via translation along a normal
    triangle at vertex 1 of the given edge embedding.

    Precondition:
    --> The given segment seg is incident to a normal triangle at vertex 1.
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


def _adjSegsAlongQuad( seg, edgeEmb, quadType, includeNonParallel ):
    """
    Yields all segments adjacent to seg via translation along the quad of the
    given type.

    If includeNonParallel is False, then only translations across parallel or
    corner faces will be taken into account. Otherwise, all translations along
    the quad will be considered.

    Precondition:
    --> The given segment seg is incident to a quad of the specified type.
    """
    tet = edgeEmb.tetrahedron()
    ver = edgeEmb.vertices()
    quadCount = seg.surface().quads( tet.index(), quadType ).pythonValue()
    triangleCount = seg.surface().triangles(
            tet.index(), ver[0] ).pythonValue()
    quadDepth = seg._segPos - triangleCount
    # From the precondition, we have 0 <= quadDepth <= quadCount.

    # The quads divide tet into two "sides". The edge opposite this
    # segment has endpoints lying on different sides, and we can label
    # these opposite endpoints opp[i], for i in {0,1}, so that ver[i] and
    # opp[i] lie on the same side of the quads, as shown in the diagram.
    #
    #                            ver[0]
    #                               •
    #                              /|\
    #                   seg.edge()/ | \
    #                            /__|__\
    #                           /|  |  |\
    #                    ver[1]•-|--|--|-•opp[1]
    #                           \|__|__|/
    #                            \  |  /
    #                             \ | /opposite edge
    #                              \|/
    #                               •
    #                            opp[0]
    #
    side = [ { 0, quadType + 1 } ]
    side.append( {0,1,2,3} - side[0] )
    if ver[0] not in side[0]:
        side[0], side[1] = side[1], side[0]
    side[0].remove( ver[0] )
    side[1].remove( ver[1] )
    opp = [ side[0].pop(), side[1].pop() ]

    # Find all edges of tet containing segments that are adjacent to seg.
    adjVertexPerms = []
    if includeNonParallel or quadDepth < quadCount:
        adjVertexPerms.append(
                Perm4( ver[0], opp[1], ver[1], opp[0] ) )
    if includeNonParallel or quadDepth > 0:
        adjVertexPerms.append(
                Perm4( opp[0], ver[1], ver[0], opp[1] ) )
    if includeNonParallel or ( quadDepth not in { 0, quadCount } ):
        adjVertexPerms.append(
                Perm4( opp[0], opp[1], ver[0], ver[1] ) )
    for vertexPerm in adjVertexPerms:
        enAdj = Edge3.faceNumber(vertexPerm)
        eiAdj = tet.edge(enAdj).index()
        verAdj = tet.edgeMapping(enAdj)
        adjTriangleCount = seg._surface.triangles(
                tet.index(), verAdj[0] ).pythonValue()
        if verAdj[0] == vertexPerm[0]:
            # Same orientation.
            yield OrientedSegment(
                    seg._surface,
                    eiAdj,
                    adjTriangleCount + quadDepth,
                    seg._orientation )
        else:
            # Opposite orientation.
            yield OrientedSegment(
                    seg._surface,
                    eiAdj,
                    adjTriangleCount + quadCount - quadDepth,
                    -1 * seg._orientation )
    return
