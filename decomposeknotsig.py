"""
Decompose a knot given by a knot signature.
"""
from sys import argv
from regina import *
from decomposeknot import decompose


if __name__ == "__main__":
    # Run decompose() with the verbose option.
    print()
    primes = decompose( Link.fromKnotSig( argv[1] ), True )
    if len(primes) == 0:
        print( "Unknot!" )
    elif len(primes) == 1:
        print( "Found 1 prime:" )
    else:
        print( "Found {} primes:".format( len(primes) ) )
    for i, loop in enumerate(primes):
        print( "    Drilled iso sig for prime #{}: {}".format(
            i, loop.drill().isoSig() ) )
