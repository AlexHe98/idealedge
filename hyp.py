"""
Use SnapPeaTriangulation and strict angle structures to attempt to certify
that an edge-ideal triangulation represents a hyperbolic 3-manifold.
"""
from regina import *
from triloops import EdgeIdealTriangulation


def knownHyperbolic(edgeIdealTri):
    """
    Is the given edge-ideal triangulation known to represent a hyperbolic
    3-manifold?

    The given edgeIdealTri should be an instance of either
    EdgeIdealTriangulation or Regina's Triangulation3.

    If this routine returns True, then the edge-ideal triangulation is
    guaranteed to represent a hyperbolic 3-manifold. On the other hand, if
    this routine returns False, we have no guarantee about whether or not the
    3-manifold is hyperbolic.

    This routine raises ValueError if the preconditions are not satisfied.

    Precondition:
    --> The triangulation is valid and orientable, and represents a
        3-manifold whose boundary is a nonempty union of tori (each of which
        is given by either an ideal loop or a real boundary torus).
    """
    if isinstance( edgeIdealTri, EdgeIdealTriangulation ):
        tri = edgeIdealTri.triangulation()
        numLoops = len(edgeIdealTri)
        drilled = edgeIdealTri.drill()
    elif isinstance( edgeIdealTri, Triangulation3 ):
        tri = edgeIdealTri
        numLoops = 0
        drilled = Triangulation3(edgeIdealTri)
    else:
        raise TypeError( "Unsupported type: {}".format(
            type(edgeIdealTri).__name__ ) )

    # Enforce the preconditions.
    if not tri.isValid():
        # This rules out:
        #   --> Invalid edges
        #   --> Vertex links which are bounded and not discs
        raise ValueError(
                "knownHyperbolic() requires a valid triangulation" )
    if not tri.isOrientable():
        raise ValueError(
                "knownHyperbolic() requires an orientable triangulation" )
    if tri.isIdeal():
        # This rules out:
        #   --> Vertex links which are closed and not 2-spheres
        raise ValueError(
                "knownHyperbolic() requires a triangulation with no ideal " +
                "vertices" )

    # At this point, we have an orientable 3-manifold such that each boundary
    # component is built entirely out of real boundary triangles. We require
    # all such boundary components to be tori.
    if numLoops + tri.countBoundaryComponents() == 0:
        raise ValueError(
                "knownHyperbolic() requires a 3-manifold with nonempty " +
                "boundary" )
    for bc in tri.boundaryComponents():
        if ( bc.eulerChar() != 0 ) or ( not bc.isOrientable() ):
            raise ValueError(
                    "knownHyperbolic() requires all real boundary " +
                    "components to be tori" )

    # If we had any ideal loops, then we already converted these to ideal
    # vertices at the beginning. Convert each real boundary torus too, and
    # then try to verify hyperbolicity.
    drilled.finiteToIdeal()
    drilled.simplify()
    simplifiedNow = True
    while simplifiedNow:
        simplifiedNow = drilled.simplify()
    spt = SnapPeaTriangulation(drilled)
    probablyHyperbolic = False
    attempts = 0
    while True:
        attempts += 1
        sol = spt.solutionType()
        try:
            # Introduced in Regina 7.4:
            geom = SnapPeaTriangulation.Solution.Geometric
            nong = SnapPeaTriangulation.Solution.Nongeometric
        except AttributeError:
            # For backwards compatibility with Regina 7.3 and earlier (but
            # this usage is deprecated as of Regina 7.4):
            geom = SnapPeaTriangulation.geometric_solution
            nong = SnapPeaTriangulation.nongeometric_solution
        if ( sol == geom or sol == nong ):
            probablyHyperbolic = True
            break
        elif attempts < 4:  # Hard-coded limit on the number of attempts.
            # Try again.
            spt.randomise()
        else:
            break
    return ( probablyHyperbolic and spt.hasStrictAngleStructure() )
