"""
Recognition of bounded orientable Seifert fibred spaces.
"""
from enum import Enum, auto
from regina import *
from aux.surface import SurfaceType, hasOnlyNonTrivialBoundaryCurves
from idealedge import ComponentDeletedByCrushing as DelComp
from idealedge import SurfaceToCrushInSuspectedSFS as CandidateSurface
from idealedge import edgeIdealTriangulationsFromCrushing
from idealedge import triangulationsWithBoundaryLoopsFromCrushing
from drill import drillMeridian
from triloops import TriangulationWithEmbeddedLoops
from triloops import EdgeIdealTriangulation, TriangulationWithBoundaryLoops


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
    ans = _recogniseVerticallyAlignedSolidTorusImpl(edgeIdealTri)
    if isinstance( ans, ManifoldProperty ):
        return ans
    fibreParams, disc = ans

    # We have a vertically-aligned solid torus if and only if crushing the
    # disc produces a triangulation which is either empty or homeomorphic to
    # the 3-ball.
    crushed = disc.crush()
    if crushed.isEmpty() or crushed.isBall():
        return fibreParams
    return ManifoldProperty.REDUCIBLE


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
        reducibility can be certified by checking that the triangulation
        given by crushing the disc D is non-empty and not homeomorphic to the
        3-ball.
    """
    if ( not edgeIdealTri.triangulation().isClosed() or
        len(edgeIdealTri) != 1 ):
        return ManifoldProperty.NOT_FST

    # Build drilled triangulation, and look for an essential disc from which
    # we can read off the fibre parameters.
    drilled = drillMeridian( edgeIdealTri[0] )
    while True:
        # Try really hard to simplify, since we will need to enumerate quad
        # vertex surfaces
        drilled.simplify()
        simplifiedNow = True
        while simplifiedNow:
            simplifiedNow = drilled.simplify()
        merEdgeIndex = drilled[0][0]

        # Search for the disc. We might find other useful surfaces instead.
        enumeration = TreeEnumeration(
                drilled.triangulation(), NormalCoords.Quad )
        while True:
            # We are enumerating finitely many surfaces, so we must
            # eventually break out of this loop.
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
                # Read off fibre parameters.
                p, q = _fibreParameters( surf, drilled )
                if p == 0:
                    return ManifoldProperty.NOT_FST
                return ( (p, q), surf )
            else:
                # We have either a 2-sphere, or a disc with trivial boundary
                # curve.
                crushAns = _crushCandidateInessentialSphereOrDisc(
                        surf, drilled )
                if crushAns == ManifoldProperty.REDUCIBLE:
                    return ManifoldProperty.REDUCIBLE

                # At this point, we should have a new drilled triangulation
                # with strictly fewer tetrahedra than before. Restart the
                # normal surface enumeration with this new triangulation.
                drilled = crushAns
                break
    raise AssertionError( "_recogniseVerticallyAlignedSolidTorusImpl() " +
                         "should never reach this point" )


def _fibreParameters( disc, drilled ):
    """
    Uses the given meridional disc to read off the Seifert fibre parameters
    carried by the given drilled triangulation.

    This routine assumes that drilled is a TriangulationWithBoundaryLoops
    whose boundary consists of exactly one real two-triangle torus, which
    contains exactly one BoundaryLoop (of length 1).
    """
    merWt = drilled.weight(disc)
    if merWt == 0:
        return (0, 1)
    merEdge = drilled.triangulation().edge( drilled[0][0] )
    front = merEdge.front()
    ver = front.vertices()
    tet = front.tetrahedron()
    lower = tet.edge( ver[0], ver[2] )
    upper = tet.edge( ver[1], ver[2] )
    lowWt = surf.edgeWeight( lower.index() ).pythonValue()
    uppWt = surf.edgeWeight( upper.index() ).pythonValue()
    if merWt == lowWt + uppWt:
        shift = lowWt
    elif uppWt == merWt + lowWt:
        shift = -lowWt
    elif lowWt == merWt + uppWt:
        shift = uppWt
    else:
        raise ValueError( "Weights don't add up." )
    q = shift % merWt
    if q > merWt // 2:
        q -= merWt
    return ( merWt, q )


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
    The specific types of the returned triangulations depends on the input:
    --> If triWithLoops is either None or an instance of
        EdgeIdealTriangulation, then each returned triangulation will be an
        instance of either EdgeIdealTriangulation or Regina's Triangulation3.
    --> If triWithLoops is an instance of TriangulationWithBoundaryLoops,
        then each returned triangulation will be an instance of either
        TriangulationWithBoundaryLoops or Regina's Triangulation3.

    Precondition
    --> The given surf should be a quadrilateral vertex normal surface.
    --> Either the given surf is a 2-sphere, or it is disc with trivial
        boundary curve.
    --> surf.triangulation() must be oriented, and each boundary component
        must be a real two-triangle torus.
    --> If triWithLoops is supplied, then triWithLoops.triangulation() should
        be the same as surf.triangulation(). In other words, triWithLoops and
        surf should both reference the same triangulation object in memory.
    """
    if ( ( triWithLoops is None ) or
        ( isinstance( triWithLoops, EdgeIdealTriangulation ) ) ):
        decomposed, _, delComps, _ = edgeIdealTriangulationsFromCrushing(
                surf, triWithLoops )
    elif isinstance( triWithLoops, TriangulationWithBoundaryLoops ):
        decomposed, delComps = triangulationsWithBoundaryLoopsFromCrushing(
                surf, triWithLoops )
    else:
        raise TypeError( "Unsupported type: {}".format(
            type(triWithLoops).__name__ ) )

    # Did we detect that the 3-manifold is reducible?
    if delComps[DelComp.L31] > 0:
        return ManifoldProperty.REDUCIBLE
    if sum( delComps.values() ) == 0 and len(decomposed) == 1:
        # The 3-manifold contains a non-separating 2-sphere.
        return ManifoldProperty.REDUCIBLE
    ans = []
    for tri in decomposed:
        # We keep only the EdgeIdealTriangulation and
        # TriangulationWithBoundaryLoop objects. Everything else should be a
        # 3-sphere, otherwise the input 3-manifold is reducible.
        if isinstance( tri, TriangulationWithEmbeddedLoops ):
            ans.append(tri)
        elif not tri.isSphere():
            return ManifoldProperty.REDUCIBLE
    return ans


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
    decomposed, numOrbCuts, delComps, inconsistent =\
            edgeIdealTriangulationsFromCrushing( surf, edgeIdealTri )

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
