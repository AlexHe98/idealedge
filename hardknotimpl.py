from overlay import overlayPD
import snappy


def attemptHardComposite(numSummands):
    # Randomly sample some knots to compose together.
    knots = snappy.HTLinkExteriors(cusps=1)
    sample = [ knots.random() for _ in range(numSummands) ]
    summands = [ k.link() for k in sample ]
    size = sum( [ len(k) for k in summands ] )
    braids = [ k.braid_word() for k in summands ]

    # Try composing.
    compPD = overlayPD(*braids)
    comp = snappy.Link(compPD)
    for _ in range(5):
        comp.simplify('global')

    # Do we get an interesting diagram?
    if ( len(comp) > 1.25*size and
        len( list( comp.dual_graph().two_cycles() ) ) == 0 ):
        #summands = comp.deconnect_sum()
        #if len(summands) != 1:
        #    continue

        # Yes! We have a "hard" diagram of a composite knot.
        newPD = comp.PD_code( min_strand_index=1 )
        return tuple( [ k.name() for k in sample ] ), newPD
    
    # No, we did not find a "hard" diagram.
    return None
