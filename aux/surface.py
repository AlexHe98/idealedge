"""
Auxiliary classes and functions for normal surfaces.
"""
from enum import Enum, auto


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
        Recognises the homeomorphism type of the given normal surface.
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
