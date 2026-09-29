"""
Decompose a knot given by an edge of a one-vertex 3-sphere.
"""
import sys
from regina import *
import snappy
from decomposeknot import decompose


if __name__ == "__main__":
    sig = sys.argv[1]
    edgeIndex = int( sys.argv[2] )
    tri = Triangulation3(sig)
    if tri.countVertices() != 1:
        print("ERROR: Triangulation must be one-vertex!")
        sys.exit()
    if not tri.isSphere():
        print("ERROR: Triangulation must be a 3-sphere!")
        sys.exit()
    primeLoops = decompose( tri.edge(edgeIndex), True )
    print()

    # Try to identify the summands.
    if hasattr( Triangulation3, "neoSig" ):
        sigGeneration = "2nd"
    else:
        #NOTE For compatibility with Regina <= 7.4.1
        sigGeneration = "1st"
    print( "Algorithm computed the following summands:" )
    for p in primeLoops:
        drilled = p.drill()
        mfd = snappy.Manifold( snappy.Triangulation( drilled.isoSig() ) )
        identified = mfd.identify()
        if identified:
            print( identified[0] )
        else:
            try:
                sig = drilled.neoSig()
            except AttributeError:
                #NOTE For compatibility with Regina <= 7.4.1
                sig = drilled.isoSig()
            print( f"Unidentified prime knot. {sigGeneration}-gen iso sig " +
                  f"for ideal triangulation: '{sig}'" )
