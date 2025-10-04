"""
Compose knot diagrams using a non-standard construction that overlays the
knots on top of each other, with the goal of obtaining a composite knot with
a diagram that is not obviously composite.
"""
from sys import argv, setrecursionlimit
import snappy
from spherogram.links.tangles import BraidTangle, join_strands
from regina import *


def compose(*knots):
    """
    Composes the given knots using a non-standard construction that overlays
    the knots on top of each other.
    """
    # The actual construction is to convert the knots into braids, and then
    # overlay the braids.
    return overlay( *[ k.braid_word() for k in knots ] )


def braidWidth(word):
    return max( abs(w) for w in word ) + 1


def overlay(*braids):
    # Create a composite knot by overlaying the given braids on each other,
    # and interleaving the strands.
    numSummands = len(braids)
    widths = [ braidWidth(b) for b in braids ]
    strands = [ (0,ii) for ii in range( widths[0] ) ]
    rightmost = [ (0,ii) for ii in range( widths[0] ) ]
    for i in range( 1, numSummands ):
        if i % 2 == 0:
            _interleaveRight( i, widths, strands, rightmost, braids[i] )
        else:
            _interleaveLeft( i, widths, strands, rightmost, braids[i] )

    # Strands have been interleaved. Now we need to introduce the crossings.
    compBraid = []
    stillProcessing = True
    while stillProcessing:
        stillProcessing = False
        for i in range( len(braids) ):
            braid = braids[i]
            if braid:
                stillProcessing = True
            else:
                continue
            oldCrossing = braid.pop(0)
            k = abs(oldCrossing)
            startStrand = strands.index( ( i, k-1 ) )
            endStrand = strands.index( ( i, k ) )
            newCrossing = ( oldCrossing // k ) * endStrand
            prefix = []
            for s in range( 1+startStrand, endStrand ):
                if strands[s][0] < i:
                    prefix.append(s)
                else:
                    prefix.append(-s)
            suffix = [ -c for c in reversed(prefix) ]
            compBraid += prefix + [newCrossing] + suffix

    # To build the desired composite knot, we take appropriate pairs of
    # strands of compBraid and join them to each other (rather than simply
    # closing up the strands like we would if we were constructing the braid
    # closure).
    pd = _overlayPD( compBraid, strands, numSummands, widths )
    return snappy.Link(pd)


def _interleaveRight( i, widths, strands, rightmost, braid ):
    # Interleave the ith braid with the previous braids that have already
    # been interleaved.
    endBand = rightmost.index( ( i-1, widths[i-1] - 1 ) )
    startBand = endBand - widths[i] + 1

    # First insert any strands that go entirely to the left of all
    # pre-existing strands of the braid.
    if startBand < 0:
        startBand *= -1
        endBand += startBand
        for ii in range(startBand):
            newStrand = ( i, ii )
            strands.insert( ii, newStrand )
            rightmost.insert( ii, newStrand )
        offset = 0
    else:
        offset = startBand

    # Now interleave the remaining strands.
    for ii in range( startBand, endBand+1 ):
        location = 1 + strands.index( rightmost[ii] )
        newStrand = ( i, ii - offset )
        strands.insert( location, newStrand )
        rightmost[ii] = newStrand

    # All done.
    return


def _interleaveLeft( i, widths, strands, rightmost, braid ):
    # Interleave the ith braid with the previous braids that have already
    # been interleaved.
    startBand = rightmost.index( ( i-1, 0 ) )
    endBand = startBand + widths[i] - 1

    # First insert strands that will actually be interleaved with
    # pre-existing strands of the braid.
    endInterleave = min( endBand+1, len(rightmost) )
    for ii in range( startBand, endInterleave ):
        location = 1 + strands.index( rightmost[ii] )
        newStrand = ( i, ii - startBand )
        strands.insert( location, newStrand )
        rightmost[ii] = newStrand

    # Now insert the remaining strands, which will go entirely to the right
    # of all pre-existing strands of the braid.
    if endBand >= endInterleave:
        for ii in range( endInterleave - startBand, widths[i] ):
            newStrand = ( i, ii )
            strands.append(newStrand)
            rightmost.append(newStrand)

    # All done.
    return


def _overlayPD( braid, threads, numSummands, widths ):
    # To build the desired composite knot, we take appropriate pairs of
    # threads of compBraid and join them to each other (rather than simply
    # closing up the threads like we would if we were constructing the braid
    # closure).
    joinedLeft = set()
    joinedRight = set()
    for i in range( numSummands - 1 ):
        if i % 2 == 0:
            j = threads.index( ( i, 0 ) )
        else:
            j = threads.index( ( i, widths[i] - 1 ) )
        joinedLeft.add(j)
        joinedRight.add(j+1)

    # Crossings are indexed in the same order as their corresponding elements
    # in the given braid word.
    totalStrands = 0
    totalCrossings = len(braid)
    pd = [ [None,None,None,None] for _ in range(totalCrossings) ]
    overcrossingSwap = set()

    # Traverse "threads" of the braid. (Here we use the word "thread" to
    # distinguish them from "strands" of the knot diagram.)
    currentThread = 0
    totalStrands += 1
    downwards = True
    while True:     # Loop to traverse threads.
        backtrack = None
        if downwards:
            # Traverse currentThread downwards.
            for i in range( len(braid) ):
                s = braid[i]

                # We have reached a crossing that exchanges threads
                # (|s| - 1) and |s|.
                if s > 0:
                    # Positive crossing.
                    if currentThread == s - 1:
                        # Undercrossing strand.
                        pd[i][0] = totalStrands
                        totalStrands += 1
                        pd[i][2] = totalStrands
                        backtrack = (i,2)
                        currentThread += 1
                        # This strand is coming in from above, so the
                        # overcrossing strand won't need to be fixed
                        # later.
                    elif currentThread == s:
                        # Overcrossing strand.
                        # We assume for now that the undercrossing strand
                        # will come in from above, and we will fix this
                        # later if necessary.
                        pd[i][3] = totalStrands
                        totalStrands += 1
                        pd[i][1] = totalStrands
                        backtrack = (i,1)
                        currentThread -= 1
                elif s < 0:
                    # Negative crossing.
                    if currentThread == -s - 1:
                        # Overcrossing strand.
                        # We assume for now that the undercrossing strand
                        # will come in from above, and we will fix this
                        # later if necessary.
                        pd[i][1] = totalStrands
                        totalStrands += 1
                        pd[i][3] = totalStrands
                        backtrack = (i,3)
                        currentThread += 1
                    elif currentThread == -s:
                        # Undercrossing strand.
                        pd[i][0] = totalStrands
                        totalStrands += 1
                        pd[i][2] = totalStrands
                        backtrack = (i,2)
                        currentThread -= 1
                        # This strand is coming in from above, so the
                        # overcrossing strand won't need to be fixed
                        # later.
                else:
                    raise ValueError()
        else:
            # Traverse currentThread upwards.
            for i in range( len(braid) - 1, -1, -1 ):
                s = braid[i]

                # We have reached a crossing that exchanges threads
                # (|s| - 1) and |s|.
                if s > 0:
                    # Positive crossing.
                    if currentThread == s - 1:
                        # Overcrossing strand.
                        # We assume for now that the undercrossing strand
                        # will come in from above, and we will fix this
                        # later if necessary.
                        pd[i][1] = totalStrands
                        totalStrands += 1
                        pd[i][3] = totalStrands
                        backtrack = (i,3)
                        currentThread += 1
                    elif currentThread == s:
                        # Undercrossing strand.
                        pd[i][0] = totalStrands
                        totalStrands += 1
                        pd[i][2] = totalStrands
                        backtrack = (i,2)
                        currentThread -= 1
                        # This strand is coming in from below, so we will
                        # need to fix the overcrossing strand later.
                        overcrossingSwap.add(i)
                elif s < 0:
                    # Negative crossing.
                    if currentThread == -s - 1:
                        # Undercrossing strand.
                        pd[i][0] = totalStrands
                        totalStrands += 1
                        pd[i][2] = totalStrands
                        backtrack = (i,2)
                        currentThread += 1
                        # This strand is coming in from below, so we will
                        # need to fix the overcrossing strand later.
                        overcrossingSwap.add(i)
                    elif currentThread == -s:
                        # Overcrossing strand.
                        # We assume for now that the undercrossing strand
                        # will come in from above, and we will fix this
                        # later if necessary.
                        pd[i][3] = totalStrands
                        totalStrands += 1
                        pd[i][1] = totalStrands
                        backtrack = (i,1)
                        currentThread -= 1
                else:
                    raise ValueError()

        # We are now at the bottom (if traversing downwards) or top (if
        # traversing upwards) of the braid. Do we turn around and join to
        # an adjacent thread, or do we continue traversing?
        if currentThread in joinedLeft:
            downwards = not downwards
            currentThread += 1
        elif currentThread in joinedRight:
            downwards = not downwards
            currentThread -= 1

        # Are we done?
        if downwards and currentThread == 0:
            # We are back to the start, so we need to backtrack and fix the
            # most recent strand.
            totalStrands -= 1
            pd[ backtrack[0] ][ backtrack[1] ] = 1

            # We might also need to fix some overcrossing strands.
            for i in overcrossingSwap:
                pd[i][1], pd[i][3] = pd[i][3], pd[i][1]

            # All done!
            return pd

    # End of traversal loop.
    # We should never reach this point.
    raise RuntimeError()


if __name__ == "__main__":
    knotNames = argv[1:]
    knots = [ snappy.Link(name) for name in knotNames ]
    composite = compose(*knots)
    print(composite)
    composite.simplify("global")
    print(composite)

    # Decompose the diagram into "diagrammatically prime" summands.
    summands = composite.deconnect_sum()
    print(summands)
    for s in summands:
        ext = s.exterior()
        print( ext.identify() )
    if len(summands) == 1:
        # We found a hard diagram of a composite knot!
        print()
        pd = composite.PD_code( min_strand_index=1 )
        print(pd)
        print()
        print( Link.fromPD(pd).knotSig() )
