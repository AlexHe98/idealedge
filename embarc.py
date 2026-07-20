"""
Embedded arcs given by sequences of segments.
"""


class EmbeddedArc:
    """
    An embedded arc given by a sequence of oriented segments.

    Typically, such arcs are constructed by splitting an embedded loop along a
    normal surface.

    Such an arc is always oriented in the sense that it will always have a
    chosen direction of traversal. The segments in the arc are indexed in
    order of traversal, and are all oriented in the direction of traversal.

    This class acts as a container of OrientedSegment objects, and therefore
    any instance embArc of this class provides the following functionality:
    --> (seg in embArc) is True if and only if seg is a segment in this arc
    --> len(embArc) is the number of segments in this arc
    --> iterating through embArc yields all the segments in order of traversal
        along the arc
    --> for i between 0 and (len(embArc) - 1), inclusive, embArc[i] returns
        the segment at index i of the arc

    This class also provides functionality to abstractly join arcs together at
    their endpoints. One use case for such joins is to indicate how arcs would
    get merged together after crushing a normal surface.
    """
    def __init__( self, segments ):
        """
        Initialises an embedded arc consisting of the given segments.

        The segments should be provided as a sequence of OrientedSegment
        objects satisfying the following conditions:
        --> The OrientedSegment objects are all defined by the same normal
            surface.
        --> The OrientedSegment objects appear in the list in the order in
            which they would be encountered as we traverse this arc.
        --> Each individual OrientedSegment is oriented in the direction of
            traversal.
        """
        self._segments = tuple(segments)
        self._joinedArc = [ None, None ]
        self._joinedEnd = [ None, None ]
        return

    def __contains__( self, segment ):
        return ( segment in self._segments )

    def __len__(self):
        return len( self._segments )

    def __iter__(self):
        return iter( self._segments )

    def __getitem__( self, key ):
        # This automatically handles both integer indices and slices.
        return self._segments.__getitem__(key)

    def endSegment( self, i ):
        """
        Returns the segment at end number i of this embedded arc.

        End 0 is the end where traversal begins, and end 1 is where traversal
        ends.

        This routine raises KeyError if i is not either 0 or 1.
        """
        if i == 0:
            return self[0]
        elif i == 1:
            return self[-1]
        raise KeyError( "Ends are numbered either 0 or 1" )

    def join( self, arcEnd, other, otherEnd ):
        """
        Abstractly joins an end of this arc to an end of the other arc.

        For both this arc and the other arc, end 0 is the end where traversal
        begins, and end 1 is where traversal ends.

        It is possible to join the two ends of this arc to each other, by
        setting other to be self, and otherEnd to be 1 - arcEnd.

        Precondition:
        --> The two ends involved in the join are not currently joined to
            anything.
        --> We are not attempting to join an edge to itself (that is, we do
            not have both other == self and otherEnd == arcEnd).

        Parameters:
        --> arcEnd      The end (either 0 or 1) of this arc to be joined to
                        the other arc.
        --> other       The other arc that will be joined to this arc at the
                        given ends.
        --> otherEnd    The end (either 0 or 1) of the other arc to which we
                        should join this arc.
        """
        self._joinImpl( arcEnd, other, otherEnd )
        other._joinImpl( otherEnd, self, arcEnd )
        return

    def _joinImpl( self, arcEnd, other, otherEnd ):
        self._joinedArc[arcEnd] = other
        self._joinedEnd[arcEnd] = otherEnd
        return

    def joinedArc( self, arcEnd ):
        """
        Returns the EmbeddedArc that is joined to this arc at the given end,
        or None if there is no such joined arc.
        """
        return self._joinedArc[arcEnd]

    def joinedEnd( self, arcEnd ):
        """
        Returns the end of self.joinedArc(arcEnd) that is joined to this arc
        at the given end, or None if there is no such join.
        """
        return self._joinedEnd[arcEnd]

    def unjoin( self, arcEnd ):
        """
        Undoes the join (if any) at the given end of this arc.

        If there was a joined arc, then this other arc will be updated
        automatically (that is, after the join has already been undone from
        this side, there is no need to call unjoin() from the other side).
        """
        other = self._joinedArc[arcEnd]
        if other is not None:
            # Before we unjoin on this side, we need to use the join data to
            # jump to the other side and unjoin there first.
            other._unjoinImpl( self._joinedEnd[arcEnd] )
            self._unjoinImpl(arcEnd)
        return

    def _unjoinImpl( self, arcEnd ):
        self._joinedArc[arcEnd] = None
        self._joinedEnd[arcEnd] = None
        return
