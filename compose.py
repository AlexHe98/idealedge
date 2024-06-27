"""
Generate random composite knots by sampling from knot tables.
"""
from sys import argv
from regina import *
from sample import sample


def generateComposite( numKnots, numSummands, *filenames ):
    """
    Generates the given number of knots, each given by composing the given
    number of summands, by randomly sampling summands from the knot table(s)
    in the given file(s).
    """
    knotNamesAndSigs = sample( numKnots * numSummands, *filenames )
    numSampled = len(knotNamesAndSigs)

    # Before generating composites, use Fisher-Yates to shuffle the sample.
    for i in range( numSampled - 1, 0, -1 ):
        swapi = RandomEngine.rand(i+1)
        knotNamesAndSigs[i], knotNamesAndSigs[swapi] = (
                knotNamesAndSigs[swapi], knotNamesAndSigs[i] )

    # Use the shuffled summands to generate a list of composite knots.
    numKnots = min( numKnots, numSampled // numSummands )
    composites = []
    for i in range(numKnots):
        compositeName, knotSig = knotNamesAndSigs[ i*numSummands ]
        compositeKnot = Link.fromKnotSig(knotSig)
        for j in range( i*numSummands + 1, (i+1)*numSummands ):
            name, knotSig = knotNamesAndSigs[j]
            compositeName += " # {}".format(name)
            compositeKnot.composeWith( Link.fromKnotSig(knotSig) )
        composites.append( ( compositeName, compositeKnot ) )
    return composites


if __name__ == "__main__":
    numKnots = int( argv[1] )
    numSummands = int( argv[2] )
    for name, _ in generateComposite( numKnots, numSummands, *argv[3:] ):
        print(name)
