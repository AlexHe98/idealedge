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
    #NOTE Both orientable() and nonOrientable() rely on the precise
    #   implementation of disc().
    if n < 1:
        raise ValueError(
                "Triangulation must have a positive number of triangles." )
    ans = Triangulation2()
    ans.newTriangles(n)
    for i in range( 1, n ):
        ans.triangle(i).join( 2, ans.triangle(i-1), Perm3(1,2) )
    return ans


def orientable( genus, boundaries ):
    """
    Returns a minimal and oriented triangulation of the orientable surface.

    Parameters:
    --> genus       The genus of the surface; this must be greater than or
                    equal to zero.
    --> boundaries  The number of boundary components in the surface; this
                    must be greater than or equal to zero.

    Returns:
        The requested orientable surface.
    """
    if genus < 0:
        raise ValueError( "Genus should be non-negative." )
    if boundaries < 0:
        raise ValueError(
                "Number of boundary components should be non-negative." )

    # For the closed case, Regina's example triangulation is already minimal.
    if boundaries == 0:
        ans = Example2.orientable( genus, 0 )

        # But it isn't guaranteed to be oriented.
        ans.orient()
        return ans

    # The disc is another special case.
    if genus == 0 and boundaries == 1:
        return disc()

    # Now construct all the other cases.
    # For the comments below, let b = boundaries and g = genus.
    if boundaries == 1:
        # We have g >= 1. The size of a minimal triangulation is 4*g - 1.
        ans = disc( 4*genus - 1 )

        #NOTE The construction below relies on the precise implementation of
        #   disc().
        for handle in range(genus):
            for faceIndex in { 4*handle, 4*handle + 1 }:
                myFace = ans.triangle(faceIndex)
                if handle == genus - 1 and faceIndex == 4*genus - 3:
                    # The very last gluing needs to be handled differently
                    # from the others.
                    yourFace = ans.triangle( 4*genus - 2 )
                    gluing = Perm3(0,1)
                else:
                    yourFace = ans.triangle( faceIndex + 2 )
                    gluing = Perm3(1,2)
                myFace.join( 0, yourFace, gluing )
    elif genus == 0:
        # We have b >= 2. The size of a minimal triangulation is 3*b - 4.
        #TODO
        raise NotImplementedError()
    else:
        # We have g >= 1 and b >= 2. The size of a minimal triangulation is
        # 4*g + 3*b - 4. We begin with the 1-boundary triangulation with
        # 4*g - 1 triangles, and then use 3*b - 3 additional triangles to
        # create all the extra boundary components.
        #TODO
        raise NotImplementedError()

    # All done!
    return ans


if __name__ == "__main__":
    availableTests = [ "disc",
                      "orbl" ]
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

    # Test orientable() routine.
    if "orbl" in testNames:
        print( "+---------------------------------+" )
        print( "| orientable( genus, boundaries ) |" )
        print( "+---------------------------------+" )
        for genus in range(9):
            #TODO Once implemented, we need to test greater number of
            #   boundary components.
            for boundaries in range(2):
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
                print()

    # If we make it here, then all tests passed.
    allTestsPassedMessage(testNames)
