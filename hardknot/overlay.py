"""
Compose knot diagrams using a non-standard construction that overlays the
knots on top of each other, with the goal of obtaining a composite knot with
a diagram that is not obviously composite.
"""
from sys import argv
from regina import *
import snappy
from hardknot.pd import PDBuilder


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
    return snappy.Link( overlayPD(*braids) )


def overlayPD( *braids, marionette=False ):
    numSummands = len(braids)
    widths = [ braidWidth(b) for b in braids ]
    strands = [ (0,ii) for ii in range( widths[0] ) ]
    compBraid = overlayBraids( *braids, marionette=marionette )

    # To build the desired composite knot, we take appropriate pairs of
    # strands of compBraid and join them to each other (rather than simply
    # closing up the strands like we would if we were constructing the braid
    # closure).
    return _overlayPDImpl( compBraid, strands, numSummands, widths )


def overlayBraids( *braids, marionette=False ):
    """
    Returns the braid constructed by overlaying the given braids and
    interleaving their strands.
    """
    # Create a composite knot by overlaying the given braids on each other,
    # and interleaving the strands.
    numSummands = len(braids)
    widths = [ braidWidth(b) for b in braids ]
    strands = [ (0,ii) for ii in range( widths[0] ) ]
    rightmost = [ (0,ii) for ii in range( widths[0] ) ]
    numCrossings = sum([ len(b) for b in braids ])
    for i in range( 1, numSummands ):
        if i % 2 == 0:
            _interleaveRight( i, widths, strands, rightmost, braids[i] )
        else:
            _interleaveLeft( i, widths, strands, rightmost, braids[i] )

    # Strands have been interleaved. Now we need to introduce the crossings.
    compBraid = []
    stillProcessing = True
    crossingsProcessed = 0
    while stillProcessing:
        stillProcessing = False
        for i in range( len(braids) ):
            braid = braids[i]
            if braid:
                stillProcessing = True
            else:
                continue

            # Marionette trick.
            if marionette:
                if crossingsProcessed == numCrossings // 4:
                    # Add in positive full twist.
                    _marionetteTwist( compBraid, len(strands), True )
                elif crossingsProcessed == 3*numCrossings // 4:
                    # Cancel out the positive full twist that we added earlier.
                    _marionetteTwist( compBraid, len(strands), False )

            # Now process crossing.
            oldCrossing = braid.pop(0)
            crossingsProcessed += 1
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
    return compBraid


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


def _marionetteTwist( braid, totalStrands, isPositive ):
    for _ in range(totalStrands):
        if isPositive:
            for s in range( totalStrands-1, 0, -1 ):
                braid.append(s)
        else:
            for s in range( 1, totalStrands ):
                braid.append(-s)
    return


def _overlayPDImpl( braid, threads, numSummands, widths ):
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
    totalCrossings = len(braid)
    builder = PDBuilder(totalCrossings)

    # Traverse "threads" of the braid. (Here we use the word "thread" to
    # distinguish them from "strands" of the knot diagram.)
    currentThread = 0
    downwards = True
    while True:     # Loop to traverse threads.
        if downwards:
            # Traverse currentThread downwards.
            for i in range(totalCrossings):
                s = braid[i]

                # We have reached a crossing that exchanges threads
                # (|s| - 1) and |s|.
                if s > 0:
                    # Positive crossing.
                    if currentThread == s - 1:
                        builder.underForwards(i)
                        currentThread += 1
                    elif currentThread == s:
                        builder.overBackwards(i)
                        currentThread -= 1
                elif s < 0:
                    # Negative crossing.
                    if currentThread == -s - 1:
                        builder.overForwards(i)
                        currentThread += 1
                    elif currentThread == -s:
                        builder.underForwards(i)
                        currentThread -= 1
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
                        builder.overForwards(i)
                        currentThread += 1
                    elif currentThread == s:
                        builder.underBackwards(i)
                        currentThread -= 1
                elif s < 0:
                    # Negative crossing.
                    if currentThread == -s - 1:
                        builder.underBackwards(i)
                        currentThread += 1
                    elif currentThread == -s:
                        builder.overBackwards(i)
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
            break
    # End of traversal loop.
    builder.finalClean()
    return builder.pd()


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
        knot = Link.fromPD(pd)
        try:
            knotSig = knot.neoSig()
        except AttributeError:
            #NOTE For compatibility with Regina <= 7.4.1
            knotSig = knot.knotSig()
            sigGeneration = "1st"
        else:
            sigGeneration = "2nd"
        print( f"{sigGeneration}-generation knot/link signature:" )
        print(knotSig)
