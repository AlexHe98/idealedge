"""
Generate a sample of hard diagrams of composite knots.
"""
from sys import argv
from timeit import default_timer
from regina import *
from hardknot import randomHardComposite


if __name__ == "__main__":
    numSamples = int( argv[1] )
    numSummands = int( argv[2] )
    workers = int( argv[3] )
    filename = argv[4]
    print( "Sampling {} hard diagrams, each with {} summands.".format(
        numSamples, numSummands ) )
    print( "{} workers.".format(workers) )
    print()

    # Save sample to file
    with ( open( filename + ".sig", "w" ) as sigs,
          open( filename + ".txt", "w" ) as pdcodes ):
        for _ in range(numSamples):
            _, pd, _ = randomHardComposite( numSummands, workers )
            pdcodes.write( str(pd) + "\n" )
            sig = Link.fromPD(pd).knotSig()
            sigs.write( sig + "\n" )
