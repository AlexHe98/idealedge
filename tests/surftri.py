"""
Test suite for minimal triangulations of surfaces.
"""
from regina import *
from surftri import orientable, nonOrientable
from sys import argv
from tests.aux import parseTestNames, allTestsPassedMessage
from tests.aux import doTest, runNamedTestSuite


def testOrientable():
    count = 0
    for genus in range(9):
        for boundaries in range(9):
            count += 1
            description = "g={}, b={}.".format( genus, boundaries )
            surf = orientable( genus, boundaries )
            if boundaries == 0:
                eulerBound = 4*genus - 2
                expectedSize = max( 2, eulerBound )
            else:
                eulerBound = 4*genus + 3*boundaries - 4
                expectedSize = max( 1, eulerBound )
            doTest( description + " Size.",
                   expectedSize, surf.size() )
            expectedEuler = 2 - 2*genus - boundaries
            doTest( description + " Euler.",
                   expectedEuler, surf.eulerChar() )
            doTest( description + " Boundary components.",
                   boundaries, surf.countBoundaryComponents() )
            doTest( description + " Oriented?", True, surf.isOriented() )
    return count


def testNonOrientable():
    count = 0
    for genus in range( 1, 9 ):
        for boundaries in range(9):
            count += 1
            description = "g={}, b={}.".format( genus, boundaries )
            surf = nonOrientable( genus, boundaries )
            if boundaries == 0:
                eulerBound = 2*genus - 2
                expectedSize = max( 2, eulerBound )
            else:
                expectedSize = 2*genus + 3*boundaries - 4
            doTest( description + " Size.",
                   expectedSize, surf.size() )
            expectedEuler = 2 - genus - boundaries
            doTest( description + " Euler.",
                   expectedEuler, surf.eulerChar() )
            doTest( description + " Boundary components.",
                   boundaries, surf.countBoundaryComponents() )
            doTest( description + " Orientable?",
                   False, surf.isOrientable() )
    return count


if __name__ == "__main__":
    availableTests = [ "orbl",
                      "norbl" ]
    testNames = parseTestNames( argv[1:], availableTests )

    # Test orientable() routine.
    if "orbl" in testNames:
        longName = "orientable( genus, boundaries )"
        runNamedTestSuite( longName, testOrientable )

    # Test nonOrientable() routine.
    if "norbl" in testNames:
        longName = "nonOrientable( genus, boundaries )"
        runNamedTestSuite( longName, testNonOrientable )

    # If we make it here, then all tests passed.
    allTestsPassedMessage(testNames)
