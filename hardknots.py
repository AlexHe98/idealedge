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
    filename = argv[3]

    # Save sample to file
    with ( open( filename + ".sig", "w" ) as sigs,
          open( filename + ".txt", "w" ) as pdcodes ):
        for _ in range(numSamples):
            _, composite = randomHardComposite(2)
            pd = composite.PD_code( min_strand_index=1 )
            pdcodes.write( str(pd) + "\n" )
            sig = Link.fromPD(pd).knotSig()
            sigs.write( sig + "\n" )
