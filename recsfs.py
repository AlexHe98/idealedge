"""
Recognition of bounded orientable Seifert fibred spaces.
"""
from enum import Enum, auto
from idealedge import ComponentDeletedByCrushing as DelComp
from idealedge import SurfaceToCrushInSuspectedSFS as CandidateSurface
from idealedge import decomposeAlong
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
