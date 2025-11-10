"""
Minimal triangulations of surfaces with nonempty boundary.
"""
from sys import argv
from regina import *
from test import parseTestNames, doTest, allTestsPassedMessage


def orientable( genus, boundaries ):
    """
    Returns a minimal and oriented triangulation of the given orientable
    surface.

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
        return _polygon(1)

    # Now construct all the other cases.
    # For the comments below, let b = boundaries and g = genus.
    if boundaries == 1:
        # We have g >= 1. The size of a minimal triangulation is 4*g - 1.
        ans = _polygon( 4*genus - 1 )

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
        ans = _polygon( 3*boundaries - 4 )

        for i in range( boundaries - 1 ):
            myFace = ans.triangle( 3*i )
            if i == boundaries - 2:
                # The very last gluing needs to be handled differently from
                # the others.
                yourFace = ans.triangle( 3*boundaries - 5 )
                gluing = Perm3(0,1)
            else:
                yourFace = ans.triangle( 3*i + 2 )
                gluing = Perm3(1,2)
            myFace.join( 0, yourFace, gluing )
    else:
        # We have g >= 1 and b >= 2. The size of a minimal triangulation is
        # 4*g + 3*b - 4. We begin with the 1-boundary triangulation with
        # 4*g - 1 triangles, and then use 3*b - 3 additional triangles to
        # create all the extra boundary components.
        ans = orientable( genus, 1 )
        _addBoundaries( ans, boundaries )

    # All done!
    return ans


def nonOrientable( genus, boundaries ):
    """
    Returns a minimal triangulation of the given non-orientable surface.

    Parameters:
    --> genus       The non-orientable genus of the surface (i.e., the number
                    of crosscaps that it contains); this must be greater than
                    or equal to one.
    --> boundaries  The number of boundary components in the surface; this
                    must be greater than or equal to zero.

    Returns:
        The requested non-orientable surface.
    """
    if genus < 1:
        raise ValueError( "Non-orientable genus should be positive." )
    if boundaries < 0:
        raise ValueError(
                "Number of boundary components should be non-negative." )

    # For the closed case, Regina's example triangulation is already minimal.
    if boundaries == 0:
        return Example2.nonOrientable( genus, 0 )

    # Now construct all the cases with nonempty boundary.
    # For the comments below, let b = boundaries and g = genus.
    if boundaries == 1:
        # We have g >= 1. The size of a minimal triangulation is 2*g - 1.
        ans = _polygon( 2*genus - 1 )

        for crosscap in range(genus):
            myFace = ans.triangle( 2*crosscap )
            if crosscap == genus - 1:
                # The very last gluing needs to be handled differently from
                # the others.
                yourFace = myFace
                gluing = Perm3(1,2,0)
            else:
                yourFace = ans.triangle( 2*crosscap + 1 )
                gluing = Perm3()
            myFace.join( 0, yourFace, gluing )
    else:
        # We have g >= 1 and b >= 2. The size of a minimal triangulation is
        # 2*g + 3*b - 4. We begin with the 1-boundary triangulation with
        # 2*g - 1 triangles, and then use 3*b - 3 additional triangles to
        # create all the extra boundary components.
        ans = nonOrientable( genus, 1 )
        _addBoundaries( ans, boundaries )

    # All done!
    return ans


def surface( genus, boundaries ):
    """
    Constructs a minimal triangulation of the surface with the given genus
    and given number of boundary components.

    The sign of the given genus parameter will be taken to indicate whether
    or not the surface is orientable. That is:
    --> for genus >= 0, the surface will be orientable; and
    --> for genus < 0, the surface will be nonorientable.
    In either case, the absolute value of the given genus parameter will be
    the genus of the surface.

    In the orientable case, the constructed triangulation will be oriented.
    """
    if genus >= 0:
        return orientable( genus, boundaries )
    else:
        return nonOrientable( -genus, boundaries )
    raise AssertionError( "This should be unreachable." )


def _polygon(n):
    """
    Returns an oriented n-triangle polygon with (n + 2) boundary edges.

    The triangulation is constructed by gluing edge (01) of triangle i to
    edge (02) of triangle (i - 1), for each i from 1 to (n - 1) (inclusive).
    """
    if n < 1:
        raise ValueError(
                "Triangulation must have a positive number of triangles." )
    ans = Triangulation2()
    ans.newTriangles(n)
    for i in range( 1, n ):
        ans.triangle(i).join( 2, ans.triangle(i-1), Perm3(1,2) )
    return ans


def _addBoundaries( surf, boundaries ):
    """
    Adds boundary components to the given one-boundary surface until it has
    the given number of boundary components.

    This routine modifies surf directly. Adding the boundary components
    increases the size of surf by (3*boundaries - 3).

    Pre-condition:
    --> surf has exactly one boundary edge, and this boundary edge is given
        by edge (01) of triangle 0.
    """
    initSize = surf.size()
    surf.insertTriangulation( _polygon( 3*boundaries - 3 ) )
    for i in range( boundaries - 1 ):
        myFace = surf.triangle( initSize + 3*i )
        yourFace = surf.triangle( initSize + 3*i + 2 )
        gluing = Perm3(1,2)
        myFace.join( 0, yourFace, gluing )
    surf.triangle(0).join( 2, surf.triangle(initSize), Perm3(0,1) )
    return


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
                print()

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
                print()

    # If we make it here, then all tests passed.
    allTestsPassedMessage(testNames)
