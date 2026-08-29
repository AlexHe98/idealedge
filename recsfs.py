"""
Recognition of bounded orientable Seifert fibred spaces.
"""
from enum import Enum, auto
from regina import *
from aux.surface import SurfaceType, hasOnlyNonTrivialBoundaryCurves
from idealedge import ComponentDeletedByCrushing as DelComp
from idealedge import SurfaceToCrushInSuspectedSFS as CandidateSurface
from idealedge import decomposeAlong
from pinch import drillMeridian
from triloops import EdgeIdealTriangulation


def recogniseSFS(tri):
    """
    Determines whether the given triangulation is a bounded orientable
    Seifert fibred space, and if so returns a Seifert fibration.

    If the triangulation is not a Seifert fibred space, then this routine
    returns None.

    This routine requires that tri is bounded and orientable, and raises
    ValueError if these conditions are not satisfied.

    Warning:
        The algorithms used in this routine rely on normal surface theory,
        and so might be very slow for larger triangulations (although faster
        tests are used where possible).
    """
    if not tri.isValid():
        # This rules out:
        #   --> Invalid edges
        #   --> Vertex links which are bounded and not discs
        raise ValueError(
                "recogniseSFS() requires a valid triangulation" )
    if not tri.isOrientable():
        raise ValueError(
                "recogniseSFS() requires an orientable triangulation" )
    if tri.isClosed() or tri.isIdeal():
        # This rules out:
        #   --> Triangulations in which all vertex links are 2-spheres
        #   --> Vertex links which are closed and not 2-spheres
        raise ValueError(
                "recogniseSFS() requires an orientable triangulation" )

    # At this point, we have an orientable 3-manifold with nonempty boundary,
    # where each boundary component is built entirely out of real boundary
    # triangles.
    for bc in tri.boundaryComponents():
        # Boundary components must all be tori to have a Seifert fibration.
        if ( bc.eulerChar() != 0 ) or ( not bc.isOrientable() ):
            return None
    tri.minimiseBoundary()

    # We now have a boundary-minimal triangulation of an orientable
    # 3-manifold whose boundary is a non-empty union of tori. This is where
    # the real work begins.

    #TODO Start with combinatorial recognition, and only fall back on normal
    #   surfaces and edge-ideal triangulations if this fails.
    #TODO Introduce a separate function which doesn't use combinatorial
    #   recognition, since this is what we want to experiment on.
    raise NotImplementedError()


class ManifoldProperty(Enum):
    """
    An enumeration of various properties that an algorithm might prove about
    an edge-ideal triangulation T of a 3-manifold M.

    The properties are not mutually exclusive, and a 3-manifold need not
    satisfy any of these properties at all.

    At present, this enumeration includes the following properties:
    --> REDUCIBLE   M is reducible.
    --> NOT_FST     T is not a vertically-aligned solid torus.
    """
    REDUCIBLE = auto()
    NOT_FST = auto()
    pass


def recogniseVerticallyAlignedSolidTorus(edgeIdealTri):
    """
    Determines whether the given EdgeIdealTriangulation is a
    vertically-aligned solid torus, and if so returns the fibre parameters
    that it carries.

    If it is not a vertically-aligned solid torus, then this routine returns
    ManifoldProperty.NOT_FST.

    Warning:
        The algorithms used in this routine rely on normal surface theory,
        and so might be very slow for larger triangulations (although faster
        tests are used where possible).
    """
    tri = edgeIdealTri.triangulation()
    if not tri.isValid() or not tri.isOrientable():
        return ManifoldProperty.NOT_FST
    #TODO
    raise NotImplementedError()


def _recogniseVerticallyAlignedSolidTorusImpl(edgeIdealTri):
    """
    Implementation of recogniseVerticallyAlignedSolidTorus().

    If this routine returns either ManifoldProperty.NOT_FST or
    ManifoldProperty.REDUCIBLE, then this is guaranteed to be correct.

    Otherwise, this routine returns candidate fibre parameters (p, q),
    together with a normal disc D which witnesses these parameters. In this
    case, we have one of the following:
    --> If the drilled 3-manifold of edgeIdealTri is irreducible, then the
        input is indeed a vertically-aligned solid torus, and moreover it
        carries a (p, q)-fibre.
    --> Otherwise, the drilled 3-manifold is not a solid torus at all, and
        boundary-reducibility can be certified by checking whether the
        triangulation given by crushing the disc D is either empty or
        homeomorphic to the 3-ball.
    """
    if ( not edgeIdealTri.triangulation().isClosed() or
        len(edgeIdealTri) != 1 ):
        return ManifoldProperty.NOT_FST

    #TODO Maybe update drillMeridian().

    # Build drilled triangulation, and look for an essential disc from which
    # we can read off the fibre parameters.
    drilled = drillMeridian( edgeIdealTri[0] )
    fibre = None    #TODO Might be enough to just use while True and break
    while fibre is None:
        # Try really hard to simplify, since we will
        # need to enumerate surfaces.
        drilled.simplify()
        simplifiedNow = True
        while simplifiedNow:
            simplifiedNow = drilled.simplify()
        merEdgeIndex = drilled[0][0]

        # Search for the disc. We might find other useful surfaces instead.
        enumeration = TreeEnumeration(
                drilled.triangulation(), NormalCoords.Quad )
        while True:
            if not enumeration.next():
                # No useful surfaces. In particular, no essential disc.
                return ManifoldProperty.NOT_FST
            surf = enumeration.buildSurface()

            # Is this a useful surface?
            surfType = SurfaceType.recognise(surf)
            if surfType == SurfaceType.RP3:
                # Orientability of the 3-manifold implies that the projective
                # plane is one-sided.
                return ManifoldProperty.REDUCIBLE
            elif surfType == SurfaceType.MOBIUS:
                if hasOnlyNonTrivialBoundaryCurves(surf):
                    # We don't consider a Mobius band with nontrivial
                    # boundary curve to be a useful surface.
                    continue
                else:
                    # Boundary of the Mobius band bounds a disc, so the
                    # 3-manifold contains a (one-sided) embedded projective
                    # plane.
                    return ManifoldProperty.REDUCIBLE
            elif surfType == SurfaceType.SPHERE:
                foundMerDisc = False
            elif surfType == SurfaceType.DISC:
                foundMerDisc = hasOnlyNonTrivialBoundaryCurves(surf)
            else:
                # Any other surface is definitely not useful.
                continue

            # Process the surface.
            if foundMerDisc:
                #TODO
                raise NotImplementedError()
            else:
                # We have either a 2-sphere, or a disc with trivial boundary
                # curve.
                #TODO
                raise NotImplementedError()
            #TODO
            raise NotImplementedError()
        #TODO
        raise NotImplementedError()
    #TODO
    raise NotImplementedError()


class _SFSpaceRecognitionInvariants:
    """
    Internal invariants used by recogniseSFS() to help recover a complete
    description of a Seifert fibration.
    """
    def __init__(self):
        """
        Default initialises the invariants.

        Specifically, we will initially have the following:
        --> self.baseEuler() will be 0.
        --> self.fibres() will be an empty list.
        --> self.isBaseNonOrientable() will be False.
        """
        self._baseEuler = 0
        self._fibres = []
        self._isBaseNonOrbl = False
        return

    def baseEuler(self):
        return self._baseEuler

    def addToBaseEuler( self, shift ):
        self._baseEuler += shift
        return

    def fibres(self):
        return self._fibres

    def newFibre( self, fibre ):
        self._fibres.append(fibre)
        return

    def isBaseNonOrientable(self):
        return self._isBaseNonOrbl

    def flagBaseNonOrientable(self):
        self._isBaseNonOrbl = True
        return


def _crushCandidateInessentialSphereOrDisc( surf, triWithLoops=None ):
    """
    Crushes the given candidate for an inessential 2-sphere or disc.

    This routine might detect that the drilled 3-manifold of the ambient
    triangulation is reducible, in which case it will return
    ManifoldProperty.REDUCIBLE.

    Otherwise, this routine returns a list consisting of the non-3-sphere
    components of the edge-ideal triangulation that results from crushing.
    Each element of the returned list will be an instance of either
    EdgeIdealTriangulation, TriangulationWithBoundaryLoops, or Regina's
    Triangulation3.
    """
    #TODO
    raise NotImplementedError()


def _crushCandidateVerticalSurface( surf, invariants, edgeIdealTri=None ):
    """
    Crushes the given candidate for a vertical surface in a
    vertically-aligned edge-ideal triangulation of a bounded orientable
    Seifert fibred space.

    This routine might detect that the drilled 3-manifold of the ambient
    triangulation is reducible, in which case it will return
    ManifoldProperty.REDUCIBLE.

    Otherwise, this routine returns a list consisting of the non-3-sphere,
    non-3-ball components of the edge-ideal triangulation that results from
    crushing. Each element of the returned list will be an instance of
    EdgeIdealTriangulation.

    The given invariants are updated in-place.

    Precondition
    --> The given surf should be a quadrilateral vertex normal surface.
    --> surf.triangulation() must be oriented.
    --> If edgeIdealTri is None, then CandidateSurface.recognise(surf) must
        be CandidateSurface.VERTICAL. Otherwise, we must have
            edgeIdealTri.weight(surf) == surf.eulerChar() > 0.
    --> If surf has real boundary, then each of its boundary curves must be a
        nontrivial curve in a two-triangle boundary torus.
    --> If edgeIdealTri is supplied, then edgeIdealTri.triangulation() should
        be the same as surf.triangulation(). In other words, edgeIdealTri and
        surf should both reference the same triangulation object in memory.
    """
    decomposed, numOrbCuts, delComps, inconsistent = decomposeAlong(
            surf, edgeIdealTri )

    # Remove all 3-sphere and 3-ball components.
    ans = []
    numSpheresAndBalls = delComps[DelComp.SPHERE] + delComps[DelComp.BALL]
    for tri in decomposed:
        if isinstance( tri, EdgeIdealTriangulation ):
            ans.append(tri)
        elif tri.isSphere() or tri.isBall():
            numSpheresAndBalls += 1
        else:
            # We have either a closed non-3-sphere, or a bounded non-3-ball.
            return ManifoldProperty.REDUCIBLE
    if numSpheresAndBalls < numOrbCuts:
        # Every cut along an orbital compression disc creates a new 2-sphere
        # boundary component, which might or might not be capped off. If the
        # input 3-manifold is irreducible, then the orbital compression discs
        # should be in bijection with 3-sphere and 3-ball components.
        return ManifoldProperty.REDUCIBLE
    assert ( numSpheresAndBalls == numOrbCuts )

    # Update the invariants.
    invariants.addToBaseEuler( numOrbCuts +
                              delComps[DelComp.FIBRE_TRIVIAL] +
                              delComps[DelComp.FIBRE_PLUS] +
                              delComps[DelComp.FIBRE_MINUS] )
    if surf.isOrientable():
        invariants.addToBaseEuler(-1)
    else:
        invariants.newFibre( (2, 1) )
    for _ in range( delComps[DelComp.FIBRE_PLUS] ):
        invariants.newFibre( (3, 1) )
    for _ in range( delComps[DelComp.FIBRE_MINUS] ):
        invariants.newFibre( (3, -1) )
    if inconsistent:
        invariants.flagBaseNonOrientable()

    # All done!
    return ans
