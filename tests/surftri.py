"""
Test suite for minimal triangulations of surfaces.
"""
from regina import *
from surftri import orientable, nonOrientable
from sys import argv
from tests.aux import parseTestNames, doTest, allTestsPassedMessage


if __name__ == "__main__":
    availableTests = [ "orbl",
                      "norbl" ]
    testNames = parseTestNames( argv[1:], availableTests )

    # Test orientable() routine.
    if "orbl" in testNames:
        print( "+---------------------------------+" )
        print( "| orientable( genus, boundaries ) |" )
        print( "+---------------------------------+" )
        for genus in range(9):
            for boundaries in range(9):
                print( "g={}, b={}".format( genus, boundaries ) )
                surf = orientable( genus, boundaries )
                if boundaries == 0:
                    eulerBound = 4*genus - 2
                    expectedSize = max( 2, eulerBound )
                else:
                    eulerBound = 4*genus + 3*boundaries - 4
                    expectedSize = max( 1, eulerBound )
                doTest( "Size.", expectedSize, surf.size() )
                expectedEuler = 2 - 2*genus - boundaries
                doTest( "Euler.", expectedEuler, surf.eulerChar() )
                doTest( "Boundary components.",
                       boundaries, surf.countBoundaryComponents() )
                doTest( "Oriented?", True, surf.isOriented() )

        # End of orientable() test.
        print()
        pass

    # Test nonOrientable() routine.
    if "norbl" in testNames:
        print( "+------------------------------------+" )
        print( "| nonOrientable( genus, boundaries ) |" )
        print( "+------------------------------------+" )
        for genus in range( 1, 9 ):
            for boundaries in range(9):
                print( "g={}, b={}".format( genus, boundaries ) )
                surf = nonOrientable( genus, boundaries )
                if boundaries == 0:
                    eulerBound = 2*genus - 2
                    expectedSize = max( 2, eulerBound )
                else:
                    expectedSize = 2*genus + 3*boundaries - 4
                doTest( "Size.", expectedSize, surf.size() )
                expectedEuler = 2 - genus - boundaries
                doTest( "Euler.", expectedEuler, surf.eulerChar() )
                doTest( "Boundary components.",
                       boundaries, surf.countBoundaryComponents() )
                doTest( "Orientable?", False, surf.isOrientable() )

        # End of nonOrientable() test.
        print()
        pass

    # If we make it here, then all tests passed.
    allTestsPassedMessage(testNames)
