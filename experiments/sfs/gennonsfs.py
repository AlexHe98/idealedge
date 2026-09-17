"""
Generate triangulations of bounded orientable 3-manifolds other than a
Seifert fibred space or hyperbolic 3-manifold.

These are intended to be representative of the worst case for the current
implementation of recogniseSFS().
"""
import sys
from timeit import default_timer
from regina import *
import snappy
from construct.sfs import orientableSFS
from recsfs import recogniseSFS, SFSRecognitionTracker


GENUS_LIST = [ -4, -3, -2, -1, 0, 1, 2 ]
FIBRE_LIST = [ (2,1), (3,1), (3,-1), (4,1), (4,-1),
              (5,1), (5,2), (5,-1), (5,-2) ]


if __name__ == "__main__":
    # Generate a random hyperbolic 3-manifold with one torus boundary
    # component, together with a random SFS with two torus boundary
    # components, and glue them together.
    tri = dict()
    hyp = snappy.OrientableCuspedCensus(num_cusps=1).random()
    tri[1] = Triangulation3( hyp.triangulation_isosig(decorated=False) )
    tri[1].idealToFinite()
    tri[1].simplify()
    tri[1].minimiseBoundary()
    RandomEngine.reseedWithHardware()
    baseSignedGenus = GENUS_LIST[ RandomEngine.rand( len(GENUS_LIST) ) ]
    numBdries = 2
    numFibres = 2 + RandomEngine.rand(3)
    fibres = [ FIBRE_LIST[ RandomEngine.rand( len(FIBRE_LIST) ) ]
              for _ in range(numFibres) ]
    tri[2] = orientableSFS( baseSignedGenus, numBdries, *fibres )
    tri[2].simplify()
    for i in [1, 2]:
#        mfd = snappy.OrientableCuspedCensus(num_cusps=i).random()
#        tri[i] = Triangulation3( mfd.triangulation_isosig(decorated=False) )
#        tri[i].idealToFinite()
#        tri[i].simplify()
#        tri[i].minimiseBoundary()
        # Try really hard to simplify.
        simplifiedNow = True
        while simplifiedNow:
            simplifiedNow = tri[i].simplify()
            if not simplifiedNow:
                simplifiedNow = tri[i].simplify()

    # For now, just arbitrarily pick an easy gluing.
    myFront = tri[1].boundaryComponent(0).edge(0).front()
    myTeti = myFront.tetrahedron().index()
    myEn = myFront.edge()
    yourFront = tri[2].boundaryComponent(0).edge(0).front()
    yourTeti = tri[1].size() + yourFront.tetrahedron().index()
    yourEn = yourFront.edge()
    tri[1].insertTriangulation( tri[2] )
    myEdge = tri[1].tetrahedron(myTeti).edge(myEn)
    yourEdge = tri[1].tetrahedron(yourTeti).edge(yourEn)
    myFront = myEdge.front()
    yourFront = yourEdge.front()
    myBack = myEdge.back()
    yourBack = yourEdge.back()
    myFront.tetrahedron().join(
            myFront.vertices()[3],
            yourFront.tetrahedron(),
            yourFront.vertices() * myFront.vertices().inverse() )
    myBack.tetrahedron().join(
            myBack.vertices()[2],
            yourBack.tetrahedron(),
            yourBack.vertices() * myBack.vertices().inverse() )

    # Might as well simplify.
    simplifiedNow = True
    while simplifiedNow:
        simplifiedNow = tri[1].simplify()
        if not simplifiedNow:
            simplifiedNow = tri[1].simplify()
    print( f"Size: {tri[1].size()}" )
    print( tri[1].isoSig() )
    print()

    # Test.
    print( "Running recogniseSFS()..." )
    print()
    sys.stdout.flush()
    useHeuristics = True
    tracker = SFSRecognitionTracker()
    start = default_timer()
    print( recogniseSFS( tri[1], useHeuristics, tracker ) )
    print( "Time: {:.6f}".format( default_timer() - start ) )
    print( "Alternate enumerations:", tracker.alternateEnumerationsCount() )
