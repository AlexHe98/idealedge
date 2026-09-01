"""
Recognition of bounded orientable Seifert fibred spaces.
"""
from enum import Enum, auto
from regina import *
from aux.looperror import BoundsDisc
from aux.surface import SurfaceType, hasOnlyNonTrivialBoundaryCurves
from hyp import knownHyperbolic
from idealedge import ComponentDeletedByCrushing as DelComp
from idealedge import SurfaceToCrushInSuspectedSFS as CandidateSurface
from idealedge import edgeIdealTriangulationsFromCrushing
from idealedge import triangulationsWithBoundaryLoopsFromCrushing
from drill import drillMeridian
from triloops import TriangulationWithEmbeddedLoops
from triloops import EdgeIdealTriangulation, TriangulationWithBoundaryLoops


def recogniseSFS( tri, useHeuristics=True ):
    """
    Determines whether the given triangulation is a bounded orientable
    Seifert fibred space, and if so returns an instance of Regina's SFSpace.

    If the triangulation is not a Seifert fibred space at all, then this
    routine returns None.

    This routine requires that tri is bounded and orientable, and raises
    ValueError if these conditions are not satisfied.

    This routine does not modify tri.

    The main Seifert fibred space recognition algorithm was designed in joint
    work with Eric Sedgwick and Jonathan Spreer; it relies on normal surface
    theory, and can therefore be very slow for larger triangulations.
    However, if useHeuristics is True (the default), this routine will
    attempt faster tests where possible. Setting useHeuristics to False is
    not recommended unless you have a particular reason for doing so, such as 
    if your goal is specifically to test the performance of the main normal
    surface algorithm.

    Warning:
        As explained above, the main algorithm used in this routine might be
        be very slow for larger triangulations.
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
                "recogniseSFS() requires a triangulation with no ideal " +
                "vertices and nonempty real boundary" )

    # At this point, we have an orientable 3-manifold with nonempty boundary,
    # where each boundary component is built entirely out of real boundary
    # triangles.
    for bc in tri.boundaryComponents():
        # Boundary components must all be tori to have a Seifert fibration.
        if ( bc.eulerChar() != 0 ) or ( not bc.isOrientable() ):
            return None
    orientedTri = Triangulation3(tri)
    orientedTri.minimiseBoundary()
    orientedTri.orient()

    # We now have an oriented boundary-minimal triangulation of a 3-manifold
    # whose boundary is a non-empty union of tori.
    #
    # Each time we come back to the top of this loop, orientedTri has
    # strictly fewer tetrahedra than before, so we can only ever loop
    # finitely many times.
    while True:
        # Try really hard to simplify, since this should help for both
        # combinatorial recognition and enumerating quad vertex surfaces.
        orientedTri.simplify()
        simplifiedNow = True
        while simplifiedNow:
            simplifiedNow = orientedTri.simplify()

        # Would like to avoid the normal surface enumeration where possible.
        if useHeuristics:
            # Attempt combinatorial recognition using Regina.
            blocked = BlockedSFS.recognise(orientedTri)
            if blocked is not None:
                return blocked.manifold()

            # Try to certify hyperbolicity.
            if knownHyperbolic(orientedTri):
                return None

            #TODO Other heuristics?

        #TODO Consider recording boundary slopes that we have already ruled
        #   out, to avoid some unnecessary computations. This might speed up
        #   cases where the input is not a Seifert fibred space.

        # Time for the heavy-duty normal surface machinery.
        #NOTE As of Regina 7.4, NS_QUAD has been deprecated, and replaced
        #   with NormalCoords.Quad.
        enumeration = TreeEnumeration( orientedTri, NormalCoords.Quad )
        while True:
            # We are enumerating finitely many surfaces, so we must
            # eventually break out of this loop.
            if not enumeration.next():
                # No vertical surfaces, so orientedTri cannot be a
                # triangulation of a Seifert fibred space.
                return None
            surf = enumeration.buildSurface()

            # Is this a useful surface?
            surfType = SurfaceType.recognise(surf)
            if surfType == SurfaceType.RP3:
                # Orientability of the 3-manifold implies that the projective
                # plane is one-sided, and hence that the 3-manifold is
                # reducible.
                return None
            elif surfType == SurfaceType.MOBIUS:
                if hasOnlyNonTrivialBoundaryCurves(surf):
                    # Mobius band with nontrivial boundary is a candidate to
                    # be a vertical surface.
                    foundCandidateVertical = True
                else:
                    # Mobius band with trivial boundary implies the existence
                    # of an embedded (one-sided) projective plane.
                    return None
            elif surfType == SurfaceType.SPHERE:
                # If the 3-manifold is a Seifert fibred space, then we expect
                # this 2-sphere to be inessential.
                foundCandidateVertical = False
            elif surfType == SurfaceType.DISC:
                if hasOnlyNonTrivialBoundaryCurves(surf):
                    # Either the 3-manifold is a solid torus, or it is
                    # reducible (and hence not a Seifert fibred space).
                    crushed = surf.crush()
                    if crushed.isEmpty() or crushed.isBall():
                        return _trivialSolidTorusFibration()
                    else:
                        # Reducible.
                        return None
                else:
                    # If the 3-manifold is a Seifert fibred space, then we
                    # expect this disc to be inessential.
                    foundCandidateVertical = False
            elif surfType == SurfaceType.ANNULUS:
                if hasOnlyNonTrivialBoundaryCurves(surf):
                    # Annulus with two nontrivial boundary curves is a
                    # candidate to be a vertical surface.
                    thin = surf.isThinEdgeLink()
                    if thin[0] is not None:
                        # Although the algorithm can handle thin edge links,
                        # in practice this usually seems unhelpful, so we
                        # just ignore them instead.
                        continue
                    foundCandidateVertical = True
                else:
                    # We don't work with annuli with trivial boundary curves.
                    continue
            else:
                # Any other surface is definitely not useful.
                continue

            # Process the surface.
            if foundCandidateVertical:
                ans = _recogniseSFSGivenCandidateVerticalSurface(surf)
                if ans is None:
                    # It turns out that the current surface is not vertical,
                    # so we need to look for another surface.
                    continue
                elif isinstance( ans, SFSpace ):
                    return ans
                elif ans == ManifoldProperty.NOT_SFS:
                    return None
                else:
                    raise AssertionError(
                            "recogniseSFS() should never reach this point" )
            else:
                crushAns = _crushCandidateInessentialSphereOrDisc(surf)
                if crushAns == ManifoldProperty.REDUCIBLE:
                    return ManifoldProperty.REDUCIBLE

                # At this point, we should have a new triangulation with
                # strictly fewer tetrahedra than before. Restart the normal
                # surface enumeration with this new triangulation.
                assert len(crushAns) == 1
                orientedTri = crushAns[0]
                assert isinstance( orientedTri, Triangulation3 )
                assert orientedTri.isOriented()
                break
        # End of enumeration loop.
    # End of loop processing triangulations.
    raise AssertionError( "recogniseSFS() should never reach this point" )


def _trivialSolidTorusFibration():
    #TODO When Regina's SFSpace is overhauled to use BundleType instead of
    #   Class, we should replace Class.bo1 with BundleType.o1.
    return SFSpace( SFSpace.Class.bo1, 0, 1 )


def recogniseTorusKnot( knot, useHeuristics=True ):
    """
    Determines whether the given knot is a torus knot, and if so returns a
    pair (p, q) describing the parameters of the torus knot.

    If the knot is not a torus knot, then this routine returns None.

    The given knot is allowed to be encoded in various ways:
    --> It could be an instance of EdgeIdealTriangulation, in which case it
        is assumed that this consists of a 3-sphere triangulation containing
        exactly one IdealLoop.
    --> It could be an instance of Regina's Edge3, in which case it is
        assumed that the endpoints of this edge are identified, and that the
        triangulation containing this edge is a 3-sphere.
    --> It could be an instance of Regina's Link or PacketOfLink, in which
        case it is assumed that this link has exactly one component.

    This routine does not modify the given knot.

    This relies on the bounded orientable Seifert fibred space recognition
    algorithm implemented in recogniseSFS(), which was designed in joint work
    with Eric Sedgwick and Jonathan Spreer. Consequently, this routine can be
    very slow for knots with many crossings, which typically require more
    tetrahedra to triangulate the exterior. If useHeuristics is True (the
    default), this routine will attempt faster tests where possible. As
    explained in the documentation of recogniseSFS(), setting useHeuristics
    to False is not recommended unless you have a particular reason for doing
    so.

    Warning:
        As explained above, the main algorithm used in this routine might be
        be very slow for knots with many crossings.
    """
    if isinstance( knot, EdgeIdealTriangulation ):
        tri = knot.drill()
    elif isinstance( knot, Edge3 ):
        tri = Triangulation3( knot.triangulation() )
        tri.pinchEdge( tri.edge( knot.index() ) )
    else:
        # Here, we assume that knot is an instance of either Link or
        # PacketOfLink.
        tri = knot.complement()
    tri.idealToFinite()
    ans = recogniseSFS( tri, useHeuristics )
    if ans is None:
        return None

    # We have a Seifert fibred space. Does it have the correct parameters to
    # be the exterior of a torus knot?
    if ans.baseClass() != SFSpace.Class.bo1:
        #TODO When Regina's SFSpace is overhauled to use BundleType instead
        #   of Class, we should replace Class.bo1 with BundleType.o1.
        return None
    if ans.baseGenus() != 0:
        return None
    if ans.punctures() != 1:
        return None
    if ans.fibreCount() != 2:
        return None

    # We have a Seifert fibred space over the disc with exactly two
    # exceptional fibres. If we have the exterior of a (p, q)-torus knot,
    # then, up to orientation, the fibre parameters can be normalised to be
    # (p, r) and (q, s) such that p*s + q*r == 1.
    myFibre = ans.fibre(0)
    yourFibre = ans.fibre(1)
    p, q, r = myFibre.alpha, yourFibre.alpha, myFibre.beta
    s = yourFibre.beta + q*ans.obstruction()
    if p > q:
        p, q = q, p
        r, s = s, r
    if (p*s + q*r) % (p*q) in { 1, p*q - 1 }:
        return (p, q)
    return None


def _recogniseSFSGivenCandidateVerticalSurface(surf):
    """
    Given a candidate vertical surface, attempts to determine whether the
    ambient triangulation is a bounded orientable Seifert fibred space.

    If recognition succeeds, then this routine returns an instance of
    Regina's SFSpace. In particular, recognition is guaranteed to succeed if
    surf is vertical with respect to *some* Seifert fibration on
    surf.triangulation().

    Otherwise, this routine returns either:
    --> ManifoldProperty.NOT_SFS, which certifies that surf.triangulation()
        does not admit any Seifert fibration; or
    --> None, which certifies that there is no Seifert fibration such that
        surf is vertical (although it is possible that some other Seifert
        fibration exists).

    Precondition
    --> The given surf is a quadrilateral vertex normal surface.
    --> CandidateSurface.recognise(surf) must be CandidateSurface.VERTICAL.
    --> surf.triangulation() is oriented, has nonempty boundary, and every
        boundary component is a real two-triangle torus.
    """
    numBdries = surf.triangulation().countBoundaryComponents()
    invariants = _SFSpaceInvariants()
    toProcess = _crushCandidateVerticalSurface( surf, invariants )
    if toProcess == ManifoldProperty.REDUCIBLE:
        return ManifoldProperty.NOT_SFS

    # At this point, toProcess is a list of EdgeIdealTriangulation objects
    # which require further processing.
    while toProcess:
        edgeIdealTri = toProcess.pop()

        # Try really hard to simplify, since we will need to enumerate quad
        # vertex surfaces
        try:
            edgeIdealTri.simplify()
            simplifiedNow = True
            while simplifiedNow:
                simplifiedNow = edgeIdealTri.simplify()
        except BoundsDisc:
            if ( len(edgeIdealTri) == 1 and
                edgeIdealTri.triangulation().isSphere() ):
                # We have found a trivial fibred solid torus.
                invariants.addToBaseEuler(1)
                continue
            else:
                # The drilled 3-manifold of edgeIdealTri is reducible.
                return ManifoldProperty.NOT_SFS

        # Search for a surface we can crush.
        #NOTE As of Regina 7.4, NS_QUAD has been deprecated, and replaced
        #   with NormalCoords.Quad.
        enumeration = TreeEnumeration(
                edgeIdealTri.triangulation(), NormalCoords.Quad )
        while True:
            # We are enumerating finitely many surfaces, so we must
            # eventually break out of this loop.
            if not enumeration.next():
                # No candidate vertical surfaces, so either edgeIdealTri is a
                # vertically-aligned solid torus, or it isn't
                # vertically-aligned at all.
                fstAns = _recogniseVerticallyAlignedSolidTorusImpl(
                        edgeIdealTri )
                if fstAns == ManifoldProperty.REDUCIBLE:
                    return ManifoldProperty.NOT_SFS
                elif fstAns == ManifoldProperty.NOT_FST:
                    return None

                # We have found a fibred solid torus.
                fibreParams, _ = fstAns
                invariants.addToBaseEuler(1)
                if fibreParams[0] > 1:
                    invariants.newFibre( SFSFibre(*fibreParams) )
                break
            surf = enumeration.buildSurface()

            # Is this a useful surface?
            surfType = SurfaceType.recognise(surf)
            wt = edgeIdealTri.weight(surf)
            if surfType == SurfaceType.RP3:
                if wt == 0:
                    # Orientability of the 3-manifold implies that the
                    # projective plane is one-sided.
                    return ManifoldProperty.NOT_SFS
                elif wt == 1:
                    # This surf restricts to a candidate vertical Mobius
                    # band.
                    foundCandidateVertical = True
                else:
                    # We do not work with higher weights.
                    continue
            elif surfType == SurfaceType.DISC:
                if hasOnlyNonTrivialBoundaryCurves(surf):
                    if wt == 0:
                        # The drilled 3-manifold is reducible.
                        return ManifoldProperty.NOT_SFS
                    elif wt == 1:
                        # This surf restricts to a candidate vertical
                        # annulus.
                        foundCandidateVertical = True
                    else:
                        # We do not work with higher weights.
                        continue
                else:
                    if wt == 0:
                        # If the drilled 3-manifold is Seifert fibred, then
                        # surf should be an inessential disc.
                        foundCandidateVertical = False
                    elif wt == 1:
                        # The drilled 3-manifold is reducible.
                        return ManifoldProperty.NOT_SFS
                    else:
                        # We do not work with higher weights.
                        continue
            elif surfType == SurfaceType.SPHERE:
                if wt == 0:
                    # If the drilled 3-manifold is Seifert fibred, then surf
                    # should be an inessential 2-sphere.
                    foundCandidateVertical = False
                elif wt == 1:
                    if ( ( len(edgeIdealTri) > 1 ) or
                        ( not edgeIdealTri.triangulation().isClosed() ) ):
                        # The drilled 3-manifold is reducible.
                        return ManifoldProperty.NOT_SFS
                    else:
                        # This edgeIdealTri is definitely not
                        # vertically-aligned.
                        return None
                elif wt == 2:
                    # This surf restricts to a candidate vertical annulus.
                    foundCandidateVertical = True
                else:
                    # We do not work with higher weights.
                    continue
            else:
                # Any other surface is definitely not useful.
                continue

            # Process the surface.
            if foundCandidateVertical:
                crushAns = _crushCandidateVerticalSurface(
                        surf, invariants, edgeIdealTri )
            else:
                crushAns = _crushCandidateInessentialSphereOrDisc(
                        surf, edgeIdealTri )
            if crushAns == ManifoldProperty.REDUCIBLE:
                return ManifoldProperty.NOT_SFS
            toProcess.extend(crushAns)
            break
        # End of enumeration loop.

    # We have emptied out toProcess, which means that the invariants carry
    # a complete description of a Seifert fibration.
    genus = 2 - invariants.baseEuler() - numBdries
    #TODO When Regina's SFSpace is overhauled to use BundleType instead of
    #   Class, we should replace Class.bo1 and Class.bn2 with BundleType.o1
    #   and BundleType.n2, respectively.
    if invariants.isBaseNonOrientable():
        baseClass = SFSpace.Class.bn2
    else:
        baseClass = SFSpace.Class.bo1
        assert ( genus % 2 == 0 )
        genus //= 2
    fibration = SFSpace( baseClass, genus, numBdries )
    for fibre in invariants.fibres():
        fibration.insertFibre(fibre)
    return fibration


class ManifoldProperty(Enum):
    """
    An enumeration of various properties that an algorithm might prove about
    an edge-ideal triangulation T of a 3-manifold M.

    The properties are not mutually exclusive, and a 3-manifold need not
    satisfy any of these properties at all.

    At present, this enumeration includes the following properties:
    --> REDUCIBLE   M is reducible.
    --> NOT_FST     T is not a vertically-aligned solid torus.
    --> NOT_SFS     T is not a vertically-aligned triangulation of a
                    Seifert fibred space
    """
    REDUCIBLE = auto()
    NOT_FST = auto()
    NOT_SFS = auto()
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
    try:
        drilled = drillMeridian( edgeIdealTri[0] )
    except BoundsDisc:
        return ManifoldProperty.NOT_FST
    while True:
        # Try really hard to simplify, since we will need to enumerate quad
        # vertex surfaces
        try:
            # Need to make sure we have minimal boundary.
            drilled.minimiseBoundary()
            simplifiedNow = True
            while simplifiedNow:
                simplifiedNow = drilled.simplify()
        except BoundsDisc:
            return ManifoldProperty.NOT_FST
        merEdgeIndex = drilled[0][0]

        # Search for the disc. We might find other useful surfaces instead.
        #NOTE As of Regina 7.4, NS_QUAD has been deprecated, and replaced
        #   with NormalCoords.Quad.
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
                assert len(crushAns) == 1
                drilled = crushAns[0]
                assert isinstance( drilled, TriangulationWithBoundaryLoops )
                assert len(drilled) == 1
                assert drilled.triangulation().countBoundaryComponents() == 1
                break
        # End of enumeration loop.
    # End of loop processing drilled triangulations.
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
    lowWt = disc.edgeWeight( lower.index() ).pythonValue()
    uppWt = disc.edgeWeight( upper.index() ).pythonValue()
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


class _SFSpaceInvariants:
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
        """
        Returns the Euler characteristic of the base surface.
        """
        return self._baseEuler

    def addToBaseEuler( self, shift ):
        """
        Adds the given (possibly negative) value to self.baseEuler().
        """
        self._baseEuler += shift
        return

    def fibres(self):
        """
        Returns the list of exceptional fibres.

        Each element of the returned list will be an instance f of Regina's
        SFSFibre with f.alpha >= 2.

        Warning:
            Modifying the returned list will invalidate the invariants which
            are being tracked by this object.
        """
        return self._fibres

    def newFibre( self, fibre ):
        """
        Adds the given new fibre to self.fibres().

        The given fibre must be an instance of Regina's SFSFibre such that
        fibre.alpha >= 2.
        """
        self._fibres.append(fibre)
        return

    def isBaseNonOrientable(self):
        """
        Returns a Boolean flag which is True if and only if the base surface
        is known to be non-orientable.
        """
        return self._isBaseNonOrbl

    def flagBaseNonOrientable(self):
        """
        Sets self.isBaseNonOrientable() to True.
        """
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
    --> The given surf is a quadrilateral vertex normal surface.
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
    --> The given surf is a quadrilateral vertex normal surface.
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
        invariants.newFibre( SFSFibre(2, 1) )
    for _ in range( delComps[DelComp.FIBRE_PLUS] ):
        invariants.newFibre( SFSFibre(3, 1) )
    for _ in range( delComps[DelComp.FIBRE_MINUS] ):
        invariants.newFibre( SFSFibre(3, -1) )
    if inconsistent:
        invariants.flagBaseNonOrientable()

    # All done!
    return ans
