"""
Try to generate a hard diagram of a composite knot by randomly composing
prime knots using the overlaying construction.
"""
from sys import setrecursionlimit
from timeit import default_timer
import snappy
from regina import *
from overlay import braidWidth, overlay


def randomHardComposite( numSummands, discriminate=True, verbose=True ):
    """
    Randomly generate a hard diagram of a composite knot by repeatedly
    composing n random prime knots, where n = numSummands.

    If discriminate is True (the default), then this routine only processes
    composite knots that are "likely" to yield a hard diagram. See the
    implementation for precise details about what this means.

    Since generating the desired hard diagram often requires many attempts
    and a significant amount of time, it is often helpful to run this routine
    in verbose mode (this is the default).
    """
    recursionErrors = 0
    if verbose:
        start = default_timer()
        attempts = 0
    while True:
        attempts += 1

        # Randomly sample some knots to compose together.
        knots = snappy.HTLinkExteriors(cusps=1)
        sample = [ knots.random() for _ in range(numSummands) ]
        summands = [ k.link() for k in sample ]
        size = sum( [ len(k) for k in summands ] )
        braids = [ k.braid_word() for k in summands ]

        # Filter out things that usually don't give a hard diagram.
        widths = [ braidWidth(b) for b in braids ]
        if discriminate:
            # Ignore cases that usually don't yield hard diagrams anyway.
            if sum(widths) >= 20:
                # This is to avoid excessive recursion due to SnapPy's
                # recursive implementation of braids and tangles.
                continue
            if max(widths) - min(widths) >= numSummands:
                # Empirically, the overlaying construction seems to be more
                # likely to lead to a hard diagram when all the overlaid
                # braids have similar widths.
                continue
        if verbose:
            print( "Attempt #{}. Time: {:.6f}. Widths: {}.".format(
                attempts, default_timer() - start, widths ) )

        # Try composing.
        try:
            comp = overlay(*braids)
        except RecursionError:
            print( "!!!!!!!!" )
            recursionErrors += 1
            if recursionErrors == 5:
                raise RuntimeError( "Too many recursion errors." )
            else:
                continue
        for _ in range(5):
            comp.simplify('global')

        # Do we get an interesting diagram?
        if ( len(comp) > 1.25*size and
            len( list( comp.dual_graph().two_cycles() ) ) == 0 ):
            #summands = comp.deconnect_sum()
            #if len(summands) != 1:
            #    continue

            # We have a "hard" diagram of a composite knot.
            if verbose:
                print()
                print( "Found a hard diagram! Time: {:.6f}.".format(
                    default_timer() - start ) )
            print()
            return [ k.name() for k in sample ], comp


if __name__ == "__main__":
    setrecursionlimit(1000000)

    # Although in principle we could try to construct hard diagrams with more
    # than 2 summands, a quirk of SnapPy's implementation of braids and
    # tangles means that we run into excessive recursion far too often to be
    # able to effectively generate hard diagrams in practice. Therefore we
    # stick to 2 summands.
    names, comp = randomHardComposite(2)
    for n in names:
        print(n)
    print()
    print(comp)
    print()
    pd = comp.PD_code( min_strand_index=1 )
    print( Link.fromPD(pd).knotSig() )
