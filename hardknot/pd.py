"""
Helper class for building PD codes.
"""


class PDBuilder:
    """
    Helper class for building a PD code.
    """
    def __init__( self, numCrossings=0 ):
        """
        Initialises a PDBuilder for building a PD code for a knot or link
        with the given number of crossings.
        """
        self._pd = [ [None,None,None,None] for _ in range(numCrossings) ]
        self._overcrossingSwap = set()
        self._totalStrands = 1
        self._backtrack = None
        return

    def newCrossing(self):
        """
        Adds a single new crossing to the PD code under construction.
        """
        self._pd.append( [None,None,None,None] )
        return

    def pd(self):
        """
        Returns the current state of the PD code under construction.

        Warning:
            You may modify the returned PD code, but such changes will modify
            the internal state of this class.
        """
        return self._pd

    def finalClean(self):
        """
        Performs a final clean-up of the constructed PD code.

        This includes fixing the last strand to be labelled 1, and fixing the
        orientation of overcrossings.

        This routine should only be called after every None value in
        self.pd() has already been replaced with a (positive) integer, and
        must only be called exactly once.
        """
        # Backtrack and fix the labelling of the last strand.
        self._totalStrands -= 1
        self._pd[ self._backtrack[0] ][ self._backtrack[1] ] = 1

        # We might also need to fix some overcrossing strands.
        for i in self._overcrossingSwap:
            self._pd[i][1], self._pd[i][3] = self._pd[i][3], self._pd[i][1]
        return

    def overForwards( self, crossingNum ):
        """
        Records the current strand as an overcrossing oriented forwards at
        the given crossing.

        This routine automatically increments the strand number by one.

        This routine assumes that the undercrossing strand will pass through
        forwards. Calling self.finalClean() at the very end of the
        construction of the PD code will automatically fix this if necessary.
        """
        self._pd[crossingNum][1] = self._totalStrands
        self._totalStrands += 1
        self._pd[crossingNum][3] = self._totalStrands
        self._backtrack = ( crossingNum, 3 )
        return

    def overBackwards( self, crossingNum ):
        """
        Records the current strand as an overcrossing oriented backwards at
        the given crossing.

        This routine automatically increments the strand number by one.

        This routine assumes that the undercrossing strand will pass through
        forwards. Calling self.finalClean() at the very end of the
        construction of the PD code will automatically fix this if necessary.
        """
        self._pd[crossingNum][3] = self._totalStrands
        self._totalStrands += 1
        self._pd[crossingNum][1] = self._totalStrands
        self._backtrack = ( crossingNum, 1 )
        return

    def underForwards( self, crossingNum ):
        """
        Records the current strand as an undercrossing oriented forwards at
        the given crossing.

        This routine automatically increments the strand number by one.
        """
        self._pd[crossingNum][0] = self._totalStrands
        self._totalStrands += 1
        self._pd[crossingNum][2] = self._totalStrands
        self._backtrack = ( crossingNum, 2 )
        return

    def underBackwards( self, crossingNum ):
        """
        Records the current strand as an undercrossing oriented backwards at
        the given crossing.

        This routine automatically increments the strand number by one.

        This routine automatically ensures that calling self.finalClean() at
        the end of the construction of the PD code will fix the orientation
        of the overcrossing.
        """
        self._pd[crossingNum][0] = self._totalStrands
        self._totalStrands += 1
        self._pd[crossingNum][2] = self._totalStrands
        self._overcrossingSwap.add(crossingNum)
        self._backtrack = ( crossingNum, 2 )
        return
