"""
Iterate through a collection of one-vertex triangulations of 3-spheres, and
perform knot decomposition experiments on all the edges of these 3-spheres.
"""
from sys import argv, stdout
from regina import *
from experiment import runDecompositionExperiment


def decomposeS3Edges(dataset):
    """
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
