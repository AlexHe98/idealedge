"""
Iterate through a collection of one-vertex triangulations of 3-spheres, and
perform knot decomposition experiments on all the edges of these 3-spheres.
"""
from sys import argv, stdout
from regina import *
from experiment import runDecompositionExperiment


def decomposeS3Edges(dataset):
    """
    Decomposes the knots represented by the edges of the one-vertex 3-spheres
    in the given dataset.

    This routine uses the readS3Edges routine to extract knots from the given
    dataset, so the format of the dataset must adhere to the specifications
    stated in the documentation for readS3Edges().

    This routine prints the results of each decomposition to standard output.
    """
    print()
    print( "+-------------------------------+" )
    print( "| Edges of one-vertex 3-spheres |" )
    print( "+-------------------------------+" )
    print()
    stdout.flush()
    slowCoefficient = 4
    runDecompositionExperiment(
            readS3Edges(dataset), slowCoefficient )
    return


def readS3Edges(dataset):
    """
    Reads edges of one-vertex 3-spheres from the given dataset.

    The given dataset should be the name of a text file with one isomorphism
    signature per line. Each listed isomorphism signature should correspond
    to a one-vertex triangulation of the 3-sphere.

    For each isomorphism signature S in the given dataset, this routine
    iterates through all the edges of the corresponding triangulation T. For
    each such edge E, this routine yields a pair (D, E), where D is a string
    that gives a description of E. Specifically, D will be of the form
    "{S}, edge {i}", where i is the index of E in T.
    """
    with open( dataset, "r" ) as lines:
        for sigline in lines:
            sig = sigline.rstrip()
            tri = Triangulation3.fromIsoSig(sig)
            for i in range( tri.countEdges() ):
                name = sig + ", edge {}".format(i)
                yield ( name, tri.edge(i) )


if __name__ == "__main__":
    decomposeS3Edges( argv[1] )
