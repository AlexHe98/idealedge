"""
Iterate through a collection of knots, given by knot signatures, and attempt
to decompose each knot into prime summands.
"""
from sys import argv, stdout
from regina import *
from experiment import runDecompositionExperiment


def decomposeKnotSigs(dataset):
    """
    Decomposes the knots represented by the knot signatures in the given
    dataset.

    This routine uses the readKnotSigs routine to extract knots from the
    given dataset, so the format of the dataset must adhere to the
    specifications stated in the documentation for readKnotSigs().

    This routine prints the results of each decomposition to standard output.
    """
    print()
    print( "+-------------------+" )
    print( "| Factorising knots |" )
    print( "+-------------------+" )
    print()
    stdout.flush()
    slowCoefficient = 2
    runDecompositionExperiment(
            readKnotSigs(dataset), slowCoefficient )
    return


def readKnotSigs(dataset):
    """
    Reads knots from the given dataset of knot signatures.

    The given dataset should be the name of a text file with one knot
    signature per line.

    For each knot signature S in the given dataset, this routine yields the
    pair (S, K), where K is the knot diagram constructed from S.
    """
    with open( dataset, "r" ) as lines:
        for sigline in lines:
            sig = sigline.rstrip()
            yield ( sig, Link.fromKnotSig(sig) )


if __name__ == "__main__":
    decomposeKnotSigs( argv[1] )
