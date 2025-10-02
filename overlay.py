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

    #TODO Remove dependence on SnapPy's recursive implementation of braids.
    # Convert compBraid into a composition of the input knots (but with a
    # more complicated diagram than the standard way to compose knots, as a
    # deliberate consequence of the above overlaying construction).
    tangle = BraidTangle(compBraid)
    top, bot = tangle.boundary
    unjoined = set( range( len(strands) ) )
    for i in range(numSummands-1):
        if i % 2 == 0:
            j = strands.index( ( i, 0 ) )
        else:
            j = strands.index( ( i, widths[i] - 1 ) )
        join_strands( tangle.adjacent[j], tangle.adjacent[j+1] )
        join_strands( tangle.adjacent[top+j], tangle.adjacent[top+j+1] )
        unjoined.remove(j)
        unjoined.remove(j+1)
    for j in unjoined:
        join_strands( tangle.adjacent[j], tangle.adjacent[top+j] )
    return snappy.Link( tangle.crossings, check_planarity=False )


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


if __name__ == "__main__":
    #TODO Remove dependence on SnapPy's recursive implementation of braids.
    # Use the overlaying construction to compose some given knots from
    # SnapPy's database, and see whether this gives a hard diagram of a
    # composite knot.
    setrecursionlimit(1000000)
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
