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

    @classmethod
    def survivors( cls, surface ):
        """
        Yields all surviving segments, with both possible orientations,
        resulting from splitting edges along the given normal surface.

        A surviving segment is an OrientedSegment that admits an associated
        surviving embedding.
        """
        #TODO
        raise NotImplementedError()

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

    def propagateTowards( self, targets ):
        #TODO Try to define "coherently"
        """
        Propagates this segment along self.surface() towards one of the given
        target segments.

        If no target segment is reachable, then this routine returns None.
        Otherwise, if the propagation reaches a target segment ts, then this
        routine returns a surviving embedding emb associated to ts such that
        emb.vertices() labels vertices 0 and 1 coherently with this segment.
        """
        #TODO
        raise NotImplementedError()

    #TODO
    pass
