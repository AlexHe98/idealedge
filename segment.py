"""
Oriented segments resulting from splitting edges along normal surfaces.
"""
from regina import *


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
            hasQuads = False
            for q in range(3):
                if surface.quads( tet.index(), q ).pythonValue() > 0:
                    hasQuads = True
                    break
            if hasQuads:
                continue

            # No quads in tet, so we will find one surviving segment along
            # each edge of tet.
            for en in range(6):
                vertexPerm = tet.edgeMapping(en)
                edgeInd = tet.edge(en).index()
                survivingSegPos = surface.triangles(
                        tet.index(), vertexPerm[0] ).pythonValue()

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

    def surface(self):
        """
        Returns the normal surface that defines this segment.
        """
        return self._surface

    def edgeIndex(self):
        """
        Returns the index of the ambient edge.
        """
        return self._edgeIndex

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
        #TODO
        raise NotImplementedError()

    def edgeWeight(self):
        """
        Returns the weight of self.surface() on the edge containing this
        segment.

        The weight will be returned as a native Python integer.
        """
        return self._surface.edgeWeight( self._edgeIndex ).pythonValue()

    def translateAlongSurface( self, targets ):
        """
        Translates this segment along self.surface(), and if possible returns
        one of the given target segments onto which this segment translates.

        If no target segment is reachable under such translation, then this
        routine returns None. Otherwise, this routine will return a target
        segment whose orientation is consistent with this segment's
        orientation under the translation.
        """
        #TODO
        raise NotImplementedError()

    #TODO
    pass
