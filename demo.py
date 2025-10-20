"""
Demo with hard diagrams of composite knots.
"""
from sys import argv
import snappy
from regina import *
from hardknot import randomHardComposite
from decomposeknot import decompose


if __name__ == "__main__":
    try:
        sig = argv[1]
    except IndexError:
        # If no knot signature is provided, then try to randomly generate a
        # hard diagram of a composite knot.
        print( "Generating a random hard diagram..." )
        print()
        numSummands = 2
        workers = 15
        summands, pd, snappyKnot = randomHardComposite( numSummands, workers )
        print("Diagram constructed from following summands:")
        for s in summands:
            print(s)
            snappy.Link(s).view()
        knot = Link.fromPD(pd)
    else:
        knot = Link.fromKnotSig(sig)
        snappyKnot = snappy.Link( knot.pdData() )

    # View the knot diagram using the PLink viewer.
    #NOTE This doesn't seem to work on all machines.
    snappyKnot.view().window.mainloop()

    # Now try to decompose the knot (in verbose mode).
    print()
    primeLoops = decompose( knot, True )
    print()

    # Finally, try to identify the summands.
    print( "Algorithm computed the following summands:" )
    last = None
    for p in primeLoops:
        drilled = p.drill()
        mfd = snappy.Manifold( snappy.Triangulation( drilled.isoSig() ) )
        identified = mfd.identify()
        if identified:
            print( identified[0] )
            last = identified[0].link().view()
        else:
            print( "Unidentified prime knot" )
    if last is not None:
        last.window.mainloop()
