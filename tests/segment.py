"""
Test suite for OrientedSegment objects.
"""
from regina import *
from segment import OrientedSegment
from sys import argv
from tests.aux import parseTestNames, allTestsPassedMessage
from tests.aux import doTest, runNamedTestSuite


def testTranslateAlongSurface():
    #TODO
    return


if __name__ == "__main__":
    availableTests = [ "surf" ]
    testNames = parseTestNames( argv[1:], availableTests )

    # Test orientable() routine.
    if "surf" in testNames:
        longName = "OrientedSegment.translateAlongSurface()"
        runNamedTestSuite( longName, testTranslateAlongSurface )

    # If we make it here, then all tests passed.
    allTestsPassedMessage(testNames)
