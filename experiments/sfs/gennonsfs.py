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


def _simplifyHard(tri):
    # Try really hard to simplify.
    simplifiedNow = True
    while simplifiedNow:
        simplifiedNow = tri.simplify()
        if not simplifiedNow:
            simplifiedNow = tri.simplify()
    return


def _randomHyp(numBdries):
    mfd = snappy.OrientableCuspedCensus(num_cusps=numBdries).random()
    tri = Triangulation3( mfd.triangulation_isosig(decorated=False) )
    tri.idealToFinite()
    tri.simplify()
    tri.minimiseBoundary()
    _simplifyHard(tri)
    return tri


def _glueTogether( myTri, yourTri, simplify=True ):
    glued = Triangulation3(myTri)
    glued.insertTriangulation( yourTri )

    # For now, just arbitrarily pick an easy gluing.
    myFront = glued.boundaryComponent(0).edge(0).front()
    myTeti = myFront.tetrahedron().index()
    myEn = myFront.edge()
    yourFront = yourTri.boundaryComponent(0).edge(0).front()
    yourTeti = myTri.size() + yourFront.tetrahedron().index()
    yourEn = yourFront.edge()
    myEdge = glued.tetrahedron(myTeti).edge(myEn)
    yourEdge = glued.tetrahedron(yourTeti).edge(yourEn)
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
    if simplify:
        _simplifyHard(glued)
    return glued


if __name__ == "__main__":
    # Build inner triangulation from random hyperbolic 3-manifold(s)
    inner = _randomHyp(1)
#    hyp = dict()
#    for hypNumBdries in [1, 2]:
#        hyp[hypNumBdries] = _randomHyp(hypNumBdries)
#    inner = _glueTogether( hyp[1], hyp[2] )

    # Glue the inner triangulation to a random SFS.
    RandomEngine.reseedWithHardware()
    baseSignedGenus = GENUS_LIST[ RandomEngine.rand( len(GENUS_LIST) ) ]
    sfsNumBdries = 2
    numFibres = 2 + RandomEngine.rand(3)
    fibres = [ FIBRE_LIST[ RandomEngine.rand( len(FIBRE_LIST) ) ]
              for _ in range(numFibres) ]
    sfsTri = orientableSFS( baseSignedGenus, sfsNumBdries, *fibres )
    _simplifyHard(sfsTri)
    tri = _glueTogether( sfsTri, inner )
    assert tri.isConnected()
    assert tri.countBoundaryComponents() == 1
    print( f"Size: {tri.size()}" )
    print( tri.isoSig() )
    print()

    # Test.
    print( "Running recogniseSFS()..." )
    print()
    sys.stdout.flush()
    useHeuristics = True
    tracker = SFSRecognitionTracker()
    start = default_timer()
    print( recogniseSFS( tri, useHeuristics, tracker ) )
    print( "Time: {:.6f}".format( default_timer() - start ) )
    print( "Alternate enumerations:", tracker.alternateEnumerationsCount() )
