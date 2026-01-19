"""
Test suite for triangulations of orientable Seifert fibre spaces.
"""
from regina import *
from surftri import surface
from sfstri import orientableSFS, OrientableBundle, TriPrism
from sys import argv, stdout
from timeit import default_timer
from tests.aux import parseTestNames, doTest, allTestsPassedMessage


if __name__ == "__main__":
    availableTests = [ "sfs",
                      "bundle",
                      "prism" ]
    testNames = parseTestNames( argv[1:], availableTests )

    # Test orientableSFS() routine.
    if "sfs" in testNames:
        print( "+-------------------------------------------------------+" )
        print( "| orientableSFS( baseSignedGenus, boundaries, *fibres ) |" )
        print( "+-------------------------------------------------------+" )
        testFibres = [
                [],
                [ (5,2) ],
                [ (3,-1), (4,1) ],
                [ (2,1), (5,1), (5,-2) ],
                [ (3,1), (7,-2), (7,2), (7,-3) ] ]
        #NOTE Regina is not always consistent with how it chooses from the
        #   many alternative names that these manifolds have.
        expectedNames = [
                "RP3 # RP3",
                "SFS [S2: (2,1) (2,1) (2,3)]",
                "SFS [RP2/n2: (3,1) (4,-1)]",
                "SFS [RP2/n2: (2,1) (5,1) (5,-2)]",
                "SFS [RP2/n2: (3,1) (7,2) (7,4) (7,-9)]",
                "M/n2 x~ S1",
                "SFS [M/n2: (5,3)]",
                "SFS [M/n2: (3,1) (4,3)]",
                "SFS [M/n2: (2,1) (5,2) (5,-1)]",
                "SFS [M/n2: (3,2) (7,2) (7,3) (7,-2)]",
                "S2 x S1",
                "RP3",
                "S3",
                "SFS [S2: (2,1) (5,1) (5,-2)]",
                "SFS [S2: (3,1) (7,2) (7,4) (7,-9)]",
                "B2 x S1",
                "B2 x S1",
                "SFS [D: (3,1) (4,-1)]",
                "SFS [D: (2,1) (5,1) (5,-2)]",
                "SFS [D: (3,1) (7,2) (7,4) (7,-9)]",
                "T x S1",
                "SFS [T: (5,2)]",
                "SFS [T: (3,1) (4,-1)]",
                "SFS [T: (2,1) (5,1) (5,-2)]",
                "SFS [T: (3,1) (7,2) (7,4) (7,-9)]",
                "Or, g=1 + 1 puncture x S1",
                "SFS [Or, g=1 + 1 puncture: (5,3)]",
                "SFS [Or, g=1 + 1 puncture: (3,1) (4,3)]",
                "SFS [Or, g=1 + 1 puncture: (2,1) (5,2) (5,-1)]",
                "SFS [Or, g=1 + 1 puncture: (3,2) (7,2) (7,3) (7,-2)]" ]

        start = default_timer()
        nameIndex = -1
        count = 0
        for genus in range( -1, 2 ):
            for boundaries in range(2):
                for fibres in testFibres:
                    description = "g={}, b={}, fibres={}.".format(
                            genus, boundaries, fibres )
                    nameIndex += 1
                    expectedName = expectedNames[nameIndex]
                    sfs = orientableSFS( genus, boundaries, *fibres )
                    doTest( description + " Oriented?",
                           True, sfs.isOriented() )
                    if expectedName == "RP3 # RP3":
                        #NOTE StandardTriangulation doesn't recognise the
                        #   minimal triangulation of RP3 # RP3.
                        print( "Skipped recognition of RP3 # RP3." )
                    else:
                        count += 1
                        actual = StandardTriangulation.recognise(sfs)
                        if actual.manifold().structure():
                            actualName = actual.manifold().structure()
                        else:
                            actualName = actual.manifold().name()
                        doTest( description + " Manifold?",
                               expectedName, actualName )

        # End of orientableSFS() test.
        print()
        print( "Tested {} orientable manifolds.".format(count) )
        print( "Time: {:.6f}".format( default_timer() - start ) )
        print( "All tests passed!" )
        print()
        stdout.flush()
        pass

    # Test OrientableBundle class.
    if "bundle" in testNames:
        print( "+------------------------+" )
        print( "| OrientableBundle class |" )
        print( "+------------------------+" )

        # Check that we get oriented triangulations of the correct
        # 3-manifolds.
        #NOTE Regina is not always consistent with how it chooses from the
        #   many alternative names that these manifolds have.
        expectedNames = [
                "KB/n2 x~ S1",
                "Non-or, g=2 + 1 puncture/n2 x~ S1",
                "Non-or, g=2 + 2 punctures/n2 x~ S1",
                "RP3 # RP3",
                "M/n2 x~ S1",
                "Non-or, g=1 + 2 punctures/n2 x~ S1",
                "S2 x S1",
                "SFS [D: (1,1)]",
                "SFS [A: (1,1)]",
                "T x S1",
                "Or, g=1 + 1 puncture x S1",
                "Or, g=1 + 2 punctures x S1",
                "Or, g=2 x S1",
                "Or, g=2 + 1 puncture x S1",
                "Or, g=2 + 2 punctures x S1" ]

        start = default_timer()
        nameIndex = -1
        count = 0
        for genus in range( -2, 3 ):
            for boundaries in range(3):
                description = "g={}, b={}.".format(
                        genus, boundaries )
                nameIndex += 1
                expectedName = expectedNames[nameIndex]
                base = surface( genus, boundaries )

                # Circle bundle.
                count += 1
                bundle = OrientableBundle( base, True )
                tri = bundle.triangulation()
                doTest( description + " Oriented?",
                       True, tri.isOriented() )
                actual = StandardTriangulation.recognise(tri).manifold()
                doTest( description + " Manifold?",
                       expectedName, actual.name() )
                for edge in base.edges():
                    if not edge.isBoundary():
                        continue
                    faceIndex = edge.front().triangle().index()
                    prism = bundle.triPrism(faceIndex)
                    square = edge.front().edge()
                    doTest(
                            description + " Square {}({}) glued?".format(
                                faceIndex, square ),
                            False, prism.isSquareGlued(square) )

        # End of OrientableBundle test.
        print()
        print( "Tested {} orientable circle bundles.".format(count) )
        print( "Time: {:.6f}".format( default_timer() - start ) )
        print( "All tests passed!" )
        print()
        stdout.flush()
        pass

    # Test TriPrism class.
    if "prism" in testNames:
        print( "+----------------+" )
        print( "| TriPrism class |" )
        print( "+----------------+" )

        def _testTri( description, tri, isSolidTorus ):
            """
            Test basic properties of the triangulation.
            """
            doTest( description + " Solid torus?",
                   isSolidTorus, tri.isSolidTorus() )
            doTest( description + " Ball?",
                   not isSolidTorus, tri.isBall() )
            doTest( description + " Oriented?",
                   True, tri.isOriented() )
            return

        def _testSlopeSign( prism, s, tri, isSolidTorus ):
            """
            Test that the slope-sign constraint holds both before and after
            flipping square s of the given prism.
            """
            if isSolidTorus:
                prismDesc = "Solid torus."
            else:
                prismDesc = "Prism."
            description = prismDesc + " Square {}.".format(s)
            # Before flipping.
            oldSlope = prism.squareSlope(s)
            slopeDesc = description + " Old slope: {}.".format(oldSlope)
            detail = " Slope-sign constraint, triangle {}."
            for t in range(2):
                doTest( slopeDesc + detail.format(t) + " Before flip.",
                       -1, oldSlope * prism.squareRoles(s,t).sign() )

            # After flipping.
            newSlope = prism.flipSlope(s)
            doTest( slopeDesc + " Flip slope.",
                    -oldSlope, newSlope )
            doTest( slopeDesc + " New slope.",
                    -oldSlope, prism.squareSlope(s) )
            _testTri( description + " After flip.", tri, isSolidTorus )
            for t in range(2):
                doTest( slopeDesc + detail.format(t) + " After flip.",
                        -1, newSlope * prism.squareRoles(s,t).sign() )
            return

        # Test both solid-torus and non-solid-torus constructions.
        start = default_timer()
        for isSolidTorus in ( True, False ):
            tri = Triangulation3()
            prism = TriPrism( tri, isSolidTorus )
            if isSolidTorus:
                description = "Solid torus. Before flip."
            else:
                description = "Prism. Before flip."
            _testTri( description, tri, isSolidTorus )

            # Flip everything one by one, and then unflip, and test that the
            # slope-sign constraint is preserved all the way through.
            for _ in range(2):
                for s in range(3):
                    _testSlopeSign( prism, s, tri, isSolidTorus )

        # End of TriPrism test.
        print()
        print( "Time: {:.6f}".format( default_timer() - start ) )
        print( "All tests passed!" )
        print()
        stdout.flush()
        pass

    # If we make it here, then all tests passed.
    allTestsPassedMessage(testNames)
