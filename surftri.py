"""
Minimal triangulations of surfaces with nonempty boundary.
"""
from sys import argv
from regina import *
from test import parseTestNames, doTest, allTestsPassedMessage


def disc(n=1):
    """
    Returns an n-triangle disc.

    By default, this returns the minimal triangulation, with n=1.

    This routine always returns an oriented triangulation of an (n+2)-sided
    polygon. In particular, this means that the returned triangulation always
    uses the minimum number of edge identifications for an n-triangle disc.
    """
    if n < 1:
        raise ValueError(
                "Triangulation must have a positive number of triangles." )
    ans = Triangulation2()
    ans.newTriangles(n)
    for i in range( 1, n ):
        ans.triangle(i).join(
                2, ans.triangle( (i-1)//2 ), Perm3( 2, i%2 ) )
    return ans


if __name__ == "__main__":
    availableTests = [ "disc" ]
    testNames = parseTestNames( argv[1:], availableTests )

    # Test disc() routine.
    if "disc" in testNames:
        print( "+---------+" )
        print( "| disc(n) |" )
        print( "+---------+" )

        # Test a range of sizes.
        for n in range( 1, 9 ):
            print( "{}-triangle disc.".format(n) )
            surf = disc(n)
            doTest( "Size.", n, surf.size() )
            doTest( "Euler.", 1, surf.eulerChar() )
            doTest( "Boundary components.",
                   1, surf.countBoundaryComponents() )
            bdryEdges = 0
            for e in surf.edges():
                if e.isBoundary():
                    bdryEdges += 1
            doTest( "Boundary edges.", n+2, bdryEdges )
            doTest( "Oriented?", True, surf.isOriented() )
            print()

    # If we make it here, then all tests passed.
    allTestsPassedMessage(testNames)
