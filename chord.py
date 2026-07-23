"""
Implements various types of embedded arcs in triangulations.
"""
from abc import ABC, abstractmethod


class NormalChord(ABC):
    """
    A "chord" for a normal surface.

    Here, a "chord" means an embedded arc whose endpoints lie on a normal
    surface S, and whose interior lies in the complement of S. This abstract
    base class does not actually implement this defining property; rather,
    such implementation is left entirely up to the concrete subclasses.

    This base class mandates the following basic functionality:
    --> A normal chord is always oriented in the sense that it will always
        have a chosen direction of traversal, and always provides a utility
        for constructing a copy with the opposite orientation.
    --> The ends of normal chords can be abstractly joined together in pairs.
        The main use case for such abstract joins is to indicate how several
        normal chords can be connected together to form abstract loops.
    """
    def __init__(self):
        """
        Initialises a normal chord with no other chords joined to its ends.
        """
        self._joinedChord = [ None, None ]
        self._joinedEnd = [ None, None ]
        return

    def join( self, chordEnd, other, otherEnd ):
        """
        Abstractly joins an end of this chord to an end of the other chord.

        For both this chord and the other chord, end 0 is the end where
        traversal begins, and end 1 is where traversal ends.

        It is possible to join the two ends of this chord to each other, by
        setting other to be self, and otherEnd to be 1 - chordEnd.

        Precondition:
        --> The two ends involved in the join are not currently joined to
            anything.
        --> We are not attempting to join an end to itself (that is, we do
            not have both other == self and otherEnd == chordEnd).

        Parameters:
        --> chordEnd    The end (either 0 or 1) of this chord to be joined to
                        the other chord.
        --> other       The other chord that will be joined to this chord at
                        the given ends.
        --> otherEnd    The end (either 0 or 1) of the other chord to which
                        we should join this chord.
        """
        self._joinImpl( chordEnd, other, otherEnd )
        other._joinImpl( otherEnd, self, chordEnd )
        return

    def _joinImpl( self, chordEnd, other, otherEnd ):
        self._joinedChord[chordEnd] = other
        self._joinedEnd[chordEnd] = otherEnd
        return

    def joinedChord( self, chordEnd ):
        """
        Returns the NormalChord that is joined to this chord at the given
        end, or None if there is no such joined chord.
        """
        return self._joinedChord[chordEnd]

    def joinedEnd( self, chordEnd ):
        """
        Returns the end of self.joinedChord(chordEnd) that is joined to this
        chord at the given end, or None if there is no such join.
        """
        return self._joinedEnd[chordEnd]

    def unjoin( self, chordEnd ):
        """
        Undoes the join (if any) at the given end of this chord.

        If there was a joined chord, then this other chord will be updated
        automatically (that is, after the join has already been undone from
        this side, there is no need to call unjoin() from the other side).
        """
        other = self._joinedChord[chordEnd]
        if other is not None:
            # Before we unjoin on this side, we need to use the join data to
            # jump to the other side and unjoin there first.
            other._unjoinImpl( self._joinedEnd[chordEnd] )
            self._unjoinImpl(chordEnd)
        return

    def _unjoinImpl( self, chordEnd ):
        self._joinedChord[chordEnd] = None
        self._joinedEnd[chordEnd] = None
        return

    @abstractmethod
    def reversed(self):
        """
        Returns a copy of this normal chord with the opposite orientation.

        The reversed copy does not retain information about how the ends are
        abstractly joined to other chords.
        """
        pass

    #TODO Other abstractmethods?
    pass


class BoundaryDualChord(NormalChord):
    """
    A normal chord which is "dual" to a boundary edge.

    This means that the chord lies entirely in the boundary of the
    triangulation, and intersects the 1-skeleton in precisely one point (the
    point of intersection therefore lies on a boundary edge which we can view
    as being dual to this chord).
    """
    def __init__( self, dualEdge, orientation ):
        """
        Initialises a boundary chord which is dual to the given boundary edge
        and which has the given orientation.

        The given dualEdge must be a boundary edge of a two-triangle torus
        boundary component of a 3-manifold triangulation.

        The given orientation must be:
        --> +1 if this chord is oriented from the front to the back of the
            embeddings of the dualEdge.
        --> -1 if this chord is oriented from the back to the front of the
            embeddings of the dualEdge.
        """
        super().__init__()  # Initialise with no other chords joined to ends.
        self._dualEdge = dualEdge
        self._orientation = orientation
        return

    def dualEdge(self):
        """
        Returns the boundary edge which is dual to this chord.
        """
        return self._dualEdge

    def orientation(self):
        """
        Returns the orientaion of this boundary dual chord.
        """
        return self._orientation

    #TODO
    pass


class SegmentedChord(NormalChord):
    """
    A normal chord built from a sequence of oriented segments.

    Typically, segmented chords are constructed by splitting an embedded loop
    along a normal surface.

    The segments in such a chord are always indexed in order of traversal,
    and are all oriented in the direction of traversal.

    In addition to providing the basic functionality of a NormalChord, this
    class acts as a container of OrientedSegment objects, and therefore any
    instance segChord of this class provides the following functionality:
    --> (seg in segChord) is True if and only if seg is a segment in this
        chord
    --> len(segChord) is the number of segments in this chord
    --> iterating through segChord yields all the segments in order of
        traversal along the chord
    --> for i between 0 and (len(segChord) - 1), inclusive, segChord[i]
        returns the segment at index i of the arc
    """
    def __init__( self, segments ):
        """
        Initialises a segmented chord built from the given sequence of
        OrientedSegment objects.

        The given segments should satisfy the following conditions:
        --> The OrientedSegment objects are all defined by the same normal
            surface.
        --> The OrientedSegment objects appear in the sequence in the order
            in which they would be encountered as we traverse this chord.
        --> Each individual OrientedSegment is oriented in the direction of
            traversal.
        """
        super().__init__()  # Initialise with no other chords joined to ends.
        self._segments = tuple(segments)
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

    def reversed(self):
        """
        Returns a copy of this segmented chord with the opposite orientation.

        The reversed copy does not retain information about how the ends are
        abstractly joined to other chords.
        """
        return SegmentedChord(
                [ seg.reversed() for seg in reversed(self._segments) ] )

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


class InternalSegmentedChord(SegmentedChord):
    """
    A segmented chord built entirely from internal segments.
    """
    #TODO
    pass


class BoundarySegmentedChord(SegmentedChord):
    """
    A segmented chord built entirely from boundary segments.
    """
    #TODO
    pass
