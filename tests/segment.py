"""
Test suite for OrientedSegment objects.
"""
from regina import *
from segment import OrientedSegment as Seg
from sys import argv
from tests.aux import parseTestNames, allTestsPassedMessage
from tests.aux import doTest, runNamedTestSuite


def testTranslateAlongSurface():
    count = 0

    # One-tetrahedron solid torus
    tri = Triangulation3("b3Na")    # 1-tet solid torus
    quadVertexSurfs = NormalSurfaces( tri, NormalCoords.Quad )
    ann, mob, disc = quadVertexSurfs
    testCases = [
            ( "Annulus, type-1 segment", Seg(ann,1,0,1),
             { Seg(ann,1,2,1), Seg(ann,1,2,-1), Seg(ann,2,0,-1), Seg(ann,2,2,1) },
             Seg(ann,1,2,-1) ),
            ( "Annulus, type-2 segment", Seg(ann,1,1,1),
             { Seg(ann,2,1,1), Seg(ann,2,1,-1) }, Seg(ann,2,1,-1) ),
            ( "Mobius, type-1 segment", Seg(mob,0,0,-1),
             { Seg(mob,0,1,-1), Seg(mob,0,1,1), Seg(mob,2,0,1), Seg(mob,2,1,-1) },
             Seg(mob,0,1,1) ),
            ( "Disc, type-1 segment", Seg(disc,0,3,1),
             { Seg(disc,1,1,-1), Seg(disc,1,1,1), Seg(disc,2,1,-1) },
             Seg(disc,1,1,1) ) ]
    for name, seg, targets, expected in testCases:
        doTest( name, expected, seg.translateAlongSurface(targets) )
        count += 1
    return count


if __name__ == "__main__":
    availableTests = [ "surf" ]
    testNames = parseTestNames( argv[1:], availableTests )

    # Test orientable() routine.
    if "surf" in testNames:
        longName = "OrientedSegment.translateAlongSurface()"
        runNamedTestSuite( longName, testTranslateAlongSurface )

    # If we make it here, then all tests passed.
    allTestsPassedMessage(testNames)
