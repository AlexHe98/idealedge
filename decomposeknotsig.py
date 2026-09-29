"""
Decompose a knot given by a knot signature.
"""
from sys import argv
from regina import *
from decomposeknot import decompose


if __name__ == "__main__":
    # Run decompose() with the verbose option.
    print()
    primes = decompose( Link.fromSig( argv[1] ), True )
    if len(primes) == 0:
        print( "Unknot!" )
    elif len(primes) == 1:
        print( "Found 1 prime:" )
    else:
        print( "Found {} primes:".format( len(primes) ) )
    if hasattr( Triangulation3, "neoSig" ):
        sigGeneration = "2nd"
    else:
        #NOTE For compatibility with Regina <= 7.4.1
        sigGeneration = "1st"
    for i, edgeIdealPrimeKnot in enumerate(primes):
        drilled = edgeIdealPrimeKnot.drill()
        try:
            sig = drilled.neoSig()
        except AttributeError:
            #NOTE For compatibility with Regina <= 7.4.1
            sig = drilled.isoSig()
        print( "    Drilled {}-gen iso sig for prime #{}: {}".format(
            sigGeneration, i, sig ) )
