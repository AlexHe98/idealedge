"""
Time randomisation of an edge-ideal triangulation.
"""
from regina import *
from loop import IdealLoop
from sys import argv, stdout
from timeit import default_timer
import cProfile


if __name__ == "__main__":
    repeats = int( argv[1] )
    sig = "7LLvLLLAvMzPLPPwPMQQAvvAMvPLzvwQMzLLQQQQcceghjmipsnouwruxyzCCyyAzxzDDDBIHIJOPNRSXYXWTU2Z5423516446635iiifclfaisnaqgiqeaaiwvcccnetaoqfonmqomumaqoofmwalmpliokscofk"
    idealEdgeIndex = 32
    tri = Triangulation3.fromIsoSig(sig)
    loop = IdealLoop( [ tri.edge(idealEdgeIndex) ] )
    for _ in range(repeats):
        start = default_timer()
        cProfile.run( "loop.randomise()" )
        print( "Time: {:.6f}.".format( default_timer() - start ) )
        print()
        stdout.flush()
