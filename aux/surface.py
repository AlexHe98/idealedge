"""
Auxiliary classes and functions for normal surfaces.
"""
from enum import Enum, auto
from regina import *


class SurfaceType(Enum):
    """
    The homeomorphism type of a normal surface.

    Currently, this enumeration recognises the following surfaces:
    --> 2-spheres (SPHERE)
    --> discs (DISC)
    --> annuli (ANNULUS)
    --> projective planes (RP3)
    --> Mobius bands (MOBIUS)
    All other surfaces will be recognised as type OTHER.
    """
    SPHERE = auto()
    DISC = auto()
    ANNULUS = auto()
    RP3 = auto()
    MOBIUS = auto()
    OTHER = auto()

    @classmethod
    def recognise( cls, surface ):
        """
        Returns the homeomorphism type of the given normal surface.
        """
        if surface.isCompact() and surface.isConnected():
            euler = surface.eulerChar()
            if surface.isOrientable():
                if euler == 2:
                    return cls.SPHERE
                if euler == 1:
                    return cls.DISC
                if euler == 0 and surface.hasRealBoundary():
                    return cls.ANNULUS
            else:
                if euler == 1:
                    return cls.RP3
                if euler == 0 and surface.hasRealBoundary():
                    return cls.MOBIUS
        # At this point, we have a surface that we don't currently recognise.
        return cls.OTHER


def hasOnlyNonTrivialBoundaryCurves(surf):
    """
    Returns True if and only if surf is either closed or has only nontrivial
    boundary curves.

    In particular, this routine vacuously returns True if surf has no
    boundary at all.

    This routine requires that surf lies inside a triangulation in which
    every boundary component is a real two-triangle torus. This requirement
    is not checked, and failure of this condition may lead to undefined
    behaviour.

    Pre-condition:
    --> Every boundary component of surf.triangulation() is a real
        two-triangle torus
    """
    tri = surf.triangulation()
    for bc in tri.boundaryComponents():
        # By assumption, bc is a real two-triangle torus. Therefore, surf has
        # a trivial boundary curve in bc if and only if every normal arc
        # coordinate in bc is positive. To check this, it suffices to examine
        # only one of the two triangles of bc.
        faceIndex = bc.triangle(0).index()
        hasTrivialBoundaryCurve = True  # True until we prove otherwise.
        for vertexNum in range(3):
            if surf.arcs( faceIndex, vertexNum ) == LargeInteger.zero:
                hasTrivialBoundaryCurve = False
                break
        if hasTrivialBoundaryCurve:
            return False
    return True


def isAnnulus(s):
    """
    Is the given normal surface s an annulus?

    Pre-condition:
    --> It is known in advance that s is connected.
    """
    return ( SurfaceType.recognise(s) == SurfaceType.ANNULUS )


def isSphere(s):
    """
    Is the given normal surface s a 2-sphere?

    Pre-condition:
    --> It is known in advance that s is connected.
    """
    return ( SurfaceType.recognise(s) == SurfaceType.SPHERE )


def countIncidentBoundaries(s):
    """
    In the triangulation containing the given normal surface s, counts the
    number of boundary components that are incident to s.

    Pre-condition:
    --> The surface s lies inside a triangulation with only real boundary
        components.
    """
    tri = s.triangulation()
    incident = set()
    for e in tri.edges():
        bdy = e.boundaryComponent()
        if ( bdy is None ) or ( bdy.index() in incident ):
            continue
        if s.edgeWeight( e.index() ).pythonValue() > 0:
            incident.add( bdy.index() )
    return len(incident)
