"""
Normal chords given by sequences of oriented segments.
"""


def pairUpChordEndsByCrushing( myChord, yourChord ):
    """
    Pairs up the ends of the two given chords according to how these ends
    would become identified after crushing the normal surface that defines
    these chords.

    This routine directly modifies the given chords by performing the
    appropriate abstract joins on their ends.

    Pre-condition:
    --> myChord and yourChord are both chords for the same normal surface.
    --> The normal surface is two-sided, and each side of the surface is
        incident to precisely two out of the four ends of the given chords.
    --> For each i in {0, 1}, both myChord.joinedChord(i) and
        yourChord.joinedChord(i) are None.
    """
    # Translate the tail of myChord along the surface. If we encounter an end
    # of yourChord, then the two chords should get joined together to form a
    # single new loop; otherwise, the two chords should be closed up to form
    # two separate loops.
    myTail = myChord.endSegment(0)
    yourTail = yourChord.endSegment(0)
    yourHead = yourChord.endSegment(1).reversed()
    targets = { yourTail, yourHead }
    targetSeg = myTail.translateAlongSurface( targets, 0 )
    if targetSeg is None:
        myChord.join( 0, myChord, 1 )
        yourChord.join( 0, yourChord, 1 )
    elif targetSeg == yourTail:
        myChord.join( 0, yourChord, 0 )
        myChord.join( 1, yourChord, 1 )
    else:   # targetSeg == yourHead
        myChord.join( 0, yourChord, 1 )
        myChord.join( 1, yourChord, 0 )
    return


class NormalChord:
    """
    A "chord" for a normal surface.

    Here, a "chord" means an embedded arc built from a sequence of oriented
    segments, such that the endpoints of the arc lie on a normal surface S,
    and the interior of the arc lies in the complement of S.

    A normal chord is always oriented in the sense that it will always have a
    chosen direction of traversal. The segments in the chord are indexed in
    order of traversal, and are all oriented in the direction of traversal.

    The ends of normal chords can be abstractly joined together in pairs. The
    main use case for such abstract joins is to indicate how several normal
    chords can be connected together to form abstract loops.

    Equality of normal chords is determined purely by their location in
    memory.

    This class also acts as a container of OrientedSegment objects. Any
    instance segChord of this class provides the following functionality:
    --> (seg in segChord) is True if and only if seg is a segment in this
        chord
    --> len(segChord) is the number of segments in this chord
    --> iterating through segChord yields all the segments in order of
        traversal along the chord
    --> for i between 0 and (len(segChord) - 1), inclusive, segChord[i]
        returns the segment at index i of the chord
    """
    def __init__( self, segments ):
        """
        Initialises a normal chord built from the given segments.

        The segments should be provided as a sequence of OrientedSegment
        objects satisfying the following conditions:
        --> The OrientedSegment objects are all defined by the same normal
            surface.
        --> The OrientedSegment objects appear in the sequence in the order
            in which they would be encountered as we traverse this chord.
        --> Each individual OrientedSegment is oriented in the direction of
            traversal.

        The chord will be initialised with no other chords abstractly joined
        to its ends.
        """
        self._segments = tuple(segments)
        self._joinedChord = [ None, None ]
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
    
    def __repr__(self):
        return "NormalChord({})".format( list(self) )

    def reversed(self):
        """
        Returns a copy of this normal chord with the opposite orientation.

        The reversed copy does not retain information about how the ends are
        abstractly joined to other chords.
        """
        return NormalChord(
                [ seg.reversed() for seg in reversed(self._segments) ] )

    def endSegment( self, i ):
        """
        Returns the segment at end number i of this normal chord.

        End 0 is the end where traversal begins, and end 1 is where traversal
        ends.

        This routine raises KeyError if i is not either 0 or 1.
        """
        if i == 0:
            return self[0]
        elif i == 1:
            return self[-1]
        raise KeyError( "Ends are numbered either 0 or 1" )

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
