"""
Experiment with certifying primeness of torus knots.
"""
from sys import argv, stdout
from math import isqrt, gcd
from regina import *
from experiment import runDecompositionExperiment


def generateTorusKnots( minCrossings, maxCrossings ):
    """
    Generates all torus knots with crossing number between minCrossings and
    maxCrossings, inclusive.

    For each generated torus knot, this routine yields a pair of the form
    (S, K), where:
    --> S is a string giving the name of the torus knot in the form
        "Torus(p, q)"; and
    --> K is a corresponding Regina Link object.
    """
    for numCrossings in range( minCrossings, maxCrossings+1 ):
        for q in range( 2, 2 + isqrt(numCrossings - 1) ):
            if numCrossings % (q-1) != 0:
                continue
            p = numCrossings // (q-1)
            if gcd(p,q) == 1:
                yield ( "Torus({}, {})".format(p,q),
                        ExampleLink.torus(p,q) )
    return


def decomposeTorusKnots( minCrossings, maxCrossings ):
    """
    Runs knot decomposition on all torus knots with crossing number between
    minCrossings and maxCrossings, inclusive.

    This routine prints the results of each decomposition to standard output.
    """
    title = "Torus knots with crossing number from {} to {}".format(
            minCrossings, maxCrossings )
    print()
    print( "+-" + "-"*len(title) + "-+" )
    print( "| {} |".format(title) )
    print( "+-" + "-"*len(title) + "-+" )
    print()
    stdout.flush()
    runDecompositionExperiment( generateTorusKnots(
        minCrossings, maxCrossings ) )
    return


if __name__ == "__main__":
    decomposeTorusKnots( int( argv[1] ), int( argv[2] ) )
