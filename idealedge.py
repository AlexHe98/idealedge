"""
Find the ideal edges after crushing a normal surface.
"""
from enum import Enum, auto
from math import gcd as pythonGCD
from regina import *
from loop import IdealLoop  #TODO Probably need to keep this, but double-check.
from triloops import EdgeIdealTriangulation
from chord import pairUpChordEndsByCrushing, NormalChord
from segment import OrientedSegment
from aux.bdry import hasOnlyMinimalRealTorusBoundaryComponents
from aux.tetrenum import tetRenumbering
from aux.quad import tetHasQuads
from aux.surface import SurfaceType, hasOnlyNonTrivialBoundaryCurves
from aux.surface import isSphere, isAnnulus, countIncidentBoundaries
from aux.looperror import NotLoop
from retriangulate.insert import layerOn
from wedge import NonSurvivingTriangularOrbitType as OrbitType
from wedge import nonSurvivingTriangularOrbitCounts as orbitCounts
#TODO Maybe move supporting classes and routines into a different file.


class SurfaceToCrushInSuspectedSFS(Enum):
    """
    An enumeration of various surfaces that can be crushed in a suspected
    real (not edge-ideal) triangulation of a bounded orientable Seifert
    fibred space (SFS).

    This enumeration includes the following types of surfaces:
    --> VERTICAL        If the ambient triangulation is a bounded orientable
                        SFS, then the surface is a candidate to be a vertical
                        annulus or Mobius band. In other words, the surface
                        is an annulus or Mobius band such that every boundary
                        curve is nontrivial.
    --> MERIDIONAL      If the ambient triangulation is a bounded orientable
                        SFS, then the surface must form a meridional disc for
                        a solid torus. In other words, the surface is a disc
                        whose boundary curve is nontrivial.
    --> BOUNDS_BALL     If the ambient triangulation is a bounded orientable
                        SFS, then the surface must bound a 3-ball. In other
                        words, either the surface is a 2-sphere, or it is a
                        disc whose boundary curve is trivial.
    """
    VERTICAL = auto()
    MERIDIONAL = auto()
    BOUNDS_BALL = auto()

    @classmethod
    def recognise( cls, surface ):
        """
        Recognises the given surface as a SurfaceToCrushInSuspectedSFS, or
        returns a string describing why the surface cannot be crushed.

        Pre-condition:
        --> Every boundary component of surface.triangulation() is a real
            two-triangle torus
        """
        surfType = SurfaceType.recognise(surface)
        if surfType == SurfaceType.ANNULUS:
            if hasOnlyNonTrivialBoundaryCurves(surface):
                return cls.VERTICAL
            return ( "To crush along an annulus, both of its boundary " +
                    "curves must be nontrivial" )
        elif surfType == SurfaceType.MOBIUS:
            if hasOnlyNonTrivialBoundaryCurves(surface):
                return cls.VERTICAL
            return ( "To crush along a Mobius band, its boundary curve " +
                    "must be nontrivial" )
        elif surfType == SurfaceType.SPHERE:
            return cls.BOUNDS_BALL
        elif surfType == SurfaceType.DISC:
            if hasOnlyNonTrivialBoundaryCurves(surface):
                return cls.MERIDIONAL
            return cls.BOUNDS_BALL
        return ( "With a real triangulation, we cannot crush along a " +
                "surface of type {}".format(surfType.name) )


class ComponentDeletedByCrushing(Enum):
    """
    An enumeration of components that are deleted by crushing because they
    arise from non-surviving triangular orbits.
    """
    BALL = auto()
    SPHERE = auto()
    L31 = auto()
    FIBRE_TRIVIAL = auto()
    FIBRE_PLUS = auto()
    FIBRE_MINUS = auto()
    pass


def decomposeAlong( surf, edgeIdealTri=None ):
    """
    Uses crushing to decompose along surf.

    This routine returns a tuple consisting of the following items:
    (0) A list of triangulated components resulting from this decomposition.
        Each element of this list will be an instance of either
        EdgeIdealTriangulation or Regina's Triangulation3. Note that this
        routine does *not* attempt to simplify the triangulations in this
        list.
    (1) A non-negative integer counting the number of orbital compression
        discs that were cut along as a consequence of crushing.
    (2) A dictionary mapping each case K of ComponentDeletedByCrushing to a
        non-negative integer counting the number of components of type K
        which arise from non-surviving triangular orbits, and which are
        therefore deleted as a consequence of crushing.
    (3) A Boolean flag which is True if and only if surf is either a
        slope-reversing annulus or a fibre-reversing 2-sphere.

    If edgeIdealTri is None, then SurfaceToCrushInSuspectedSFS.recognise(surf)
    must return an instance of SurfaceToCrushInSuspectedSFS. Otherwise,
    edgeIdealTri.allowsCrush(surf) must be True. This routine raises
    ValueError if surf does not satisfy these conditions.

    We also require surf to be a quadrilateral vertex normal surface, but
    this routine does *not* check this condition.

    Precondition
    --> The given surf should be a quadrilateral vertex normal surface.
    --> If surf has real boundary, then each boundary component that it meets
        must be a two-triangle torus.
    --> If edgeIdealTri is supplied, then edgeIdealTri.triangulation() should
        be the same as surf.triangulation(). In other words, edgeIdealTri and
        surf should both reference the same triangulation object in memory.
    """
    tri = surf.triangulation()
    if not hasOnlyMinimalRealTorusBoundaryComponents(tri):
        raise ValueError( "decomposeAlong() requires that the ambient " +
                         "triangulation only has real boundary components " +
                         "that are two-triangle tori" )

    # Compute the sequences of chords which will become new ideal loops (if
    # any) after crushing surf. Along the way, we also check that we are
    # actually allowed to crush surf.
    if edgeIdealTri is None:
        crushCase = SurfaceToCrushInSuspectedSFS.recognise(surf)
        if isinstance( crushCase, str ):
            # Not allowed to crush.
            raise ValueError(crushCase)
        elif crushCase == SurfaceToCrushInSuspectedSFS.VERTICAL:
            # We have boundary chords which piece together to form new ideal
            # loops after crushing.
            chordSequences = _buildNewLoopsFromBoundaryChords(surf)
        else:
            # We are crushing a 2-sphere or disc, so we don't pick up any new
            # ideal loops after crushing.
            chordSequences = []
    else:
        # Enforce the precondition that the two input objects reference
        # precisely the same Triangulation3 object in memory.
        if tri is not edgeIdealTri.triangulation():
            raise RuntimeError( "decomposeAlong() requires the input " +
                               "NormalSurface and the input " +
                               "EdgeIdealTriangulation to reference the " +
                               "same Triangulation3 object in memory" )
        edgeIdealTri.checkCrushAllowed(surf)
        chordSequences = _buildNewLoopsFromIdealChords( surf, edgeIdealTri )

    # Convert chord sequences into sequences of surviving edge embeddings.
    newLoops = []
    numOrbitalCompressions = 0
    numNonSurvivingLoops = 0
    foundInconsistentLoop = False
    survivors = OrientedSegment.survivors(surf)
    for chordsInNewLoop, loopStatus in chordSequences:
        if loopStatus == _IdealLoopStatus.COMPRESSED:
            numOrbitalCompressions += 1
            # Loop is compressed away, so no need to add it to the newLoops.
            continue
        elif loopStatus == _IdealLoopStatus.INCONSISTENT:
            foundInconsistentLoop = True

        # As a consequence of the pre-condition that surf is quad vertex,
        # together with the fact this new loop is not compressed away by an
        # orbital compression disc, we know that every segment of this loop
        # must be simple. It follows that one segment of this loop survives
        # if and only if every segment survives. In other words, the
        # pre-condition for _extractSurvivingEmbeddings() is satisfied.
        survivingEmbs = _extractSurvivingEmbeddings(
                chordsInNewLoop, survivors )
        if survivingEmbs:
            newLoops.append(survivingEmbs)
        else:
            numNonSurvivingLoops += 1

    # Count deleted components arising from non-surviving triangular orbits.
    deletedOrbitCounts = orbitCounts(surf)
    DelComp = ComponentDeletedByCrushing
    deletedComponentCounts = {
            DelComp.BALL: deletedOrbitCounts[OrbitType.BOUNDARY] }
    if numNonSurvivingLoops:
        # From the pre-condition that surf is quad vertex, every wedge cycle
        # with a twist corresponds to a deleted fibre of multiplicity 3.
        deletedComponentCounts[ DelComp.L31 ] = 0
        twistPlus = deletedOrbitCounts[ OrbitType.TWIST_PLUS ]
        twistMinus = deletedOrbitCounts[ OrbitType.TWIST_MINUS ]
        deletedComponentCounts[ DelComp.FIBRE_PLUS ] = twistPlus
        deletedComponentCounts[ DelComp.FIBRE_MINUS ] = twistMinus

        # The number of deleted 3-spheres, vs deleted fibres of multiplicity
        # 1, is then determined by numNonSurvivingLoops.
        trivFibres = numNonSurvivingLoops - twistPlus - twistMinus
        assert trivFibres >= 0
        deletedComponentCounts[ DelComp.FIBRE_TRIVIAL ] = trivFibres
        delSpheres = ( deletedOrbitCounts[
            OrbitType.TRIVIAL_CYCLE ] - trivFibres )
        assert delSpheres >= 0
        deletedComponentCounts[ DelComp.SPHERE ] = delSpheres
    else:
        # No deleted fibres.
        deletedComponentCounts[ DelComp.FIBRE_PLUS ] = 0
        deletedComponentCounts[ DelComp.FIBRE_MINUS ] = 0
        deletedComponentCounts[ DelComp.FIBRE_TRIVIAL ] = 0
        deletedComponentCounts[ DelComp.L31 ] =\
                ( deletedOrbitCounts[ OrbitType.TWIST_PLUS ] +
                 deletedOrbitCounts[ OrbitType.TWIST_MINUS ] )
        deletedComponentCounts[ DelComp.SPHERE ] =\
                ( deletedOrbitCounts[ OrbitType.TRIVIAL_CYCLE ] )

    # Find where the new ideal loops will be after first crushing, and then
    # closing up pinched boundary 2-spheres.
    doomed = [ tet for tet in surf.triangulation().tetrahedra()
              if tetHasQuads( surf, tet.index() ) ]
    tetIndicesAfterCrush = tetRenumbering(doomed)
    crushed = surf.crush()
    loopEmbSeqs = []
    for oldEmbSequence in newLoops:
        embSequence = []
        i = 0
        while i < len(oldEmbSequence):
            oldEmb = oldEmbSequence[i]
            crushedTet = crushed.tetrahedron(
                    tetIndicesAfterCrush[ oldEmb.tetrahedron().index() ] )
            crushedEdge = crushedTet.edge( oldEmb.edge() )
            crushedEdgeMapping = crushedTet.edgeMapping( oldEmb.edge() )
            tail = crushedEdgeMapping.inverse()[ oldEmb.vertices()[0] ]
            head = 1 - tail
            if crushedEdge.isBoundary():
                # Because of the promises made by _findBoundaryChords(), we
                # should have two consecutive boundary edges. After layering
                # (if necessary) and closing up, we replace these two edges
                # with a single edge.
                i += 1
                assert i < len(oldEmbSequence)
                front = crushedEdge.front()
                back = crushedEdge.back()
                if ( front.tetrahedron().triangle( front.vertices()[3] ) ==
                    back.tetrahedron().triangle( back.vertices()[2] ) ):
                    # Need to do a layering before closing up.
                    layerEdge = front.tetrahedron().edge(
                            front.vertices()[head], front.vertices()[2] )
                    layerTet = layerOn(layerEdge)
                    closedEdge = layerTet.edge(5)
                    closedHead = layerTet.edgeMapping(5).inverse()[
                            front.tetrahedron().adjacentFace(
                                front.vertices()[3] ) ]
                else:
                    # No layering needed before closing up.
                    closedEdgeNum = Edge3.faceNumber(
                            front.vertices()[tail], front.vertices()[2] )
                    closedEdge = front.tetrahedron().edge(closedEdgeNum)
                    closedHead = front.tetrahedron().edgeMapping(
                            closedEdgeNum ).inverse()[ front.vertices()[2] ]

                # Perform the closing up.
                closedFront = closedEdge.front()
                closedBack = closedEdge.back()
                closedFront.tetrahedron().join(
                        closedFront.vertices()[3],
                        closedBack.tetrahedron(),
                        closedBack.vertices() * Perm4(2, 3) *
                        closedFront.vertices().inverse() )
                if closedHead == 0:
                    newEmb = EdgeEmbedding3(
                            closedFront.tetrahedron(),
                            closedFront.vertices() * Perm4(1, 0, 3, 2) )
                else:
                    newEmb = closedFront
            else:
                newEmb = EdgeEmbedding3( crushedTet, oldEmb.vertices() )
            embSequence.append(newEmb)
            i += 1
        loopEmbSeqs.append(embSequence)

    # Split crushed into its components.
    if crushed.isConnected():
        components = [crushed]
        compLoopInfo = [[]]
        for embSequence in loopEmbSeqs:
            compLoopInfo[0].append(
                    [ ( emb.tetrahedron().index(), emb.vertices() )
                     for emb in embSequence ] )
    else:
        components = list( crushed.triangulateComponents() )

        # Work out how tetrahedra get renumbered after splitting crushed into
        # its components.
        shiftedIndex = []
        compSize = [0] * crushed.countComponents()
        for i in range( crushed.size() ):
            compi = crushed.tetrahedron(i).component().index()
            shiftedIndex.append( compSize[compi] )
            compSize[compi] += 1

        # Using the renumbering that we just computed, record shifted
        # tetrahedron indices for the ideal loops.
        compLoopInfo = [ [] for _ in range( crushed.countComponents() ) ]
        for embSequence in loopEmbSeqs:
            singleLoopInfo = []
            for survivingEmb in embSequence:
                teti = survivingEmb.tetrahedron().index()
                singleLoopInfo.append( ( shiftedIndex[teti],
                                        survivingEmb.vertices() ) )

            # Abuse the fact that teti persists beyond the scope of the above
            # for loop.
            compi = crushed.tetrahedron(teti).component().index()
            compLoopInfo[compi].append(singleLoopInfo)

    # Use compLoopInfo to find the ideal loops in each component.
    triList = []
    for compi in range( crushed.countComponents() ):
        tri = components[compi]
        loopInfo = compLoopInfo[compi]
        loops = []
        for singleLoopInfo in loopInfo:
            # To construct an IdealLoop, we need:
            #   --> a list of edges, in order as we traverse the loop; and
            #   --> an orientation, which is either +1 if the first edge of
            #       the loop is oriented from vertex 0 to vertex 1, and -1 if
            #       the first edge is oriented from vertex 1 to vertex 0.
            edgeList = []
            for teti, ver in singleLoopInfo:
                edgeList.append(
                        tri.tetrahedron(teti).edge( ver[0], ver[1] ) )
            firstTet = tri.tetrahedron( singleLoopInfo[0][0] )
            firstVer = singleLoopInfo[0][1]
            firstEdgeNum = Edge3.faceNumber(firstVer)
            if firstVer[0] == firstTet.edgeMapping(firstEdgeNum)[0]:
                orientation = 1
            else:
                orientation = -1
            loops.append( IdealLoop( edgeList, orientation ) )

        # If we have any loops at all, then package them all together as a
        # single EdgeIdealTriangulation. Otherwise, just add an ordinary
        # Triangulation3 to the triList.
        if loops:
            triList.append( EdgeIdealTriangulation(loops) )
        else:
            triList.append(tri)

    # All done!
    return ( triList, numOrbitalCompressions, deletedComponentCounts,
            foundInconsistentLoop )


class _IdealLoopStatus(Enum):
    """
    Status of a new ideal loop created by crushing a quad vertex surface.

    This is one of the following:
    COMPRESSED      Indicates that the loop consists entirely of two type-1
                    segments which cobound an orbital compression disc, which
                    means that the loop gets compressed away as a side-effect
                    of crushing
    CONSISTENT      Indicates that the loop is not compressed, and that the
                    normal chords comprising the loop are all consistently
                    oriented
    INCONSISTENT    Indicates that the loop is not compressed, but includes
                    some pair of normal chords that are inconsistently
                    oriented
    """
    COMPRESSED = auto()
    CONSISTENT = auto()
    INCONSISTENT = auto()
    pass


def _buildNewLoopsFromBoundaryChords(surf):
    """
    Uses boundary chords to build the new ideal loops that would arise from
    crushing the given normal surface surf.

    In detail, this routine returns a list, each of whose elements describes
    a new ideal loop via a pair consisting of the following items:
    (0) A list of boundary chords, appearing in order of traversal around the
        new loop, and also oriented consistently with the order of traversal.
    (1) A status given by _IdealLoopStatus.

    Warning:
        This routine does not check any of the pre-conditions listed below.

    Pre-condition:
    --> The given surf should be a quadrilateral vertex normal surface.
    --> The ambient triangulation surf.triangulation() must have nonempty
        minimal toroidal boundary.
    --> SurfaceToCrushInSuspectedSFS.recognise(surf) must be
        SurfaceToCrushInSuspectedSFS.VERTICAL.
    """
    tri = surf.triangulation()
    boundaryChords = _findBoundaryChords(surf)

    # Put the boundary chords together to form the new loops. From the
    # pre-conditions, the total number of boundary chords is equal to the
    # number of boundary curves of surf, which is either 1 or 2.
    chordSequences = []
    if len(boundaryChords) == 2:
        myChord, yourChord = boundaryChords
        pairUpChordEndsByCrushing( myChord, yourChord )
        if myChord.joinedChord(0) == yourChord:
            # The two boundary chords join together to form a single loop.
            # This means that surf is an annulus which either is
            # slope-reversing or spans two distinct boundary tori. From the
            # pre-condition that surf is quad-vertex, together with the
            # promises made by _findBoundaryChords(), we know that:
            #   (a) the loop is not compressed away by an orbital compression
            #       disc; and
            #   (b) surf is slope-reversing if and only if one of the two
            #       boundary chords has length 1.
            chordsInNewLoop = [myChord]
            chordOrientationsDisagree = ( myChord.joinedEnd(0) == 0 )
            if chordOrientationsDisagree:
                chordsInNewLoop.append( yourChord.reversed() )
            else:
                chordsInNewLoop.append(yourChord)
            if ( len(myChord) == 1 or len(yourChord) == 1 ):
                assert chordOrientationsDisagree
                loopStatus = _IdealLoopStatus.INCONSISTENT
            else:
                # In this case, we are joining two boundary chords from two
                # different boundary tori, so the choices of orientations on
                # the two chords are arbitrary and meaningless.
                loopStatus = _IdealLoopStatus.CONSISTENT
            chordSequences.append( ( chordsInNewLoop, loopStatus ) )
            return chordSequences

    # At this point, we know that each boundary chord simply joins with
    # itself to form its own loop.
    #
    # In this case, it is necessary to check for orbital compression discs.
    # Because _findBoundaryChords() returns length-2 boundary chords whenever
    # possible, we know that if an orbital compression disc D exists, then we
    # will have a boundary chord which bounds D.
    for bdryChord in boundaryChords:
        if _boundsOrbitalCompressionDisc(bdryChord):
            loopStatus = _IdealLoopStatus.COMPRESSED
        else:
            loopStatus = _IdealLoopStatus.CONSISTENT
        chordSequences.append( ( [bdryChord], loopStatus ) )
    return chordSequences


def _buildNewLoopsFromIdealChords( surf, edgeIdealTri ):
    """
    Uses the ideal chords to build the new ideal loops that would arise from
    crushing the given normal surface surf.

    In detail, this routine returns a list, each of whose elements describes
    a new ideal loop via a pair consisting of the following items:
    (0) A list of normal chords, appearing in order of traversal around the
        new loop, and also oriented consistently with the order of traversal.
    (1) A status given by _IdealLoopStatus.

    The new ideal loops described by the returned list are related to the old
    ideal loops in edgeIdealTri as follows:
    --> Old ideal loops that are disjoint from the surface will be left
        topologically untouched. In particular, their orientations will be
        preserved.
    --> Old ideal loops that intersect the surface will be split into
        multiple ideal chords; additionally, if the surface is a disc, then
        there might be a single boundary chord. Each such (ideal or boundary)
        chord either survives the crushing operation, or is entirely
        destroyed by crushing. For the surviving chords, crushing will
        essentially rearrange how the endpoints of these chords are joined
        together, thereby yielding new ideal loops.

    Warning:
        This routine does not check any of the pre-conditions listed below.

    Pre-condition:
    --> The given surf should be a quadrilateral vertex normal surface.
    --> The ambient triangulation surf.triangulation() must either be closed
        or have minimal toroidal boundary.
    --> Both surf.triangulation() and edgeIdealTri.triangulation() must
        reference the same Triangulation3 object in memory.
    --> edgeIdealTri.allowsCrush(surf) must be True.
    """
    tri = surf.triangulation()
    idealChords = edgeIdealTri.splitIntoChords(surf)
    boundaryChords = _findBoundaryChords(surf)

    # Put the ideal chords together to form the new loops. If one of the
    # ideal chords has unjoined ends, then these ends will need to be joined
    # to the ends of the (unique) boundary chord.
    chordSequences = []
    while idealChords:
        currentChord = idealChords.pop()
        chordsInNewLoop = [currentChord]
        loopStatus = _IdealLoopStatus.CONSISTENT    # Until proven otherwise.
        currentTailEnd = currentChord.joinedEnd(1)
        currentChord = currentChord.joinedChord(1)
        if currentChord is None:
            # From the pre-conditions, surf must be a disc with edge-ideal
            # weight 1, and we must have found the unique ideal chord which
            # needs to be joined with the unique boundary chord.
            assert len(boundaryChords) == 1
            bdryChord = boundaryChords.pop()
            pairUpChordEndsByCrushing( chordsInNewLoop[0], bdryChord )

            # The current choice of orientation on the boundary chord is
            # arbitrary and meaningless. We re-orient (if necessary) to
            # ensure that the orientation is consistent with the pre-existing
            # ideal chord.
            if chordsInNewLoop[0].joinedEnd(0) == 0:
                chordsInNewLoop.append( bdryChord.reversed() )
            else:
                chordsInNewLoop.append(bdryChord)
        else:
            # Traverse the new loop, and pick up all its constituent ideal
            # chords.
            while currentChord != chordsInNewLoop[0]:
                idealChords.remove(currentChord)
                if currentTailEnd == 1:
                    # The current chord is oriented inconsistently with the
                    # first chord in this new loop.
                    loopStatus = _IdealLoopStatus.INCONSISTENT
                    chordsInNewLoop.append( currentChord.reversed() )
                else:
                    chordsInNewLoop.append(currentChord)
                currentHeadEnd = 1 - currentTailEnd

                # Move on to the next chord in the loop.
                currentTailEnd = currentChord.joinedEnd(currentHeadEnd)
                currentChord = currentChord.joinedChord(currentHeadEnd)
            assert currentTailEnd == 0  # Tail of the first chord

        # Check whether this loop is compressed away by an orbital
        # compression disc.
        if len(chordsInNewLoop) == 1:
            assert loopStatus == _IdealLoopStatus.CONSISTENT
            if _boundsOrbitalCompressionDisc( chordsInNewLoop[0] ):
                loopStatus = _IdealLoopStatus.COMPRESSED

        # Done with this loop.
        chordSequences.append( ( chordsInNewLoop, loopStatus ) )

    # All done!
    return chordSequences


def _extractSurvivingEmbeddings( chordSequence, survivors ):
    """
    Translates the sequence of normal chords along parallel cells to obtain a
    sequence of surviving edge embeddings.

    This routine returns a (possibly empty) list of edge embeddings, where
    each such edge embedding is given by calling seg.survivingEmbedding(),
    for some OrientedSegment seg in the given set of survivors.

    This routine will never modify the given chordSequence and the given set
    of survivors.

    Pre-condition:
    --> If some segment of some chord in chordSequence translates along
        parallel cells to a surviving segment, then the same holds for all
        such segments.
    """
    newLoopEmbeddings = []
    for chord in chordSequence:
        for seg in chord:
            # Note that OrientedSegment.translateAlongParallelCells()
            # promises to never modify survivors.
            survivingSeg = seg.translateAlongParallelCells(survivors)
            if survivingSeg is None:
                # Using the pre-condition, we may conclude that the entire
                # loop doesn't survive.
                return []
            newLoopEmbeddings.append(
                    survivingSeg.survivingEmbedding() )
    return newLoopEmbeddings


def _findBoundaryChords(surf):
    """
    Returns a set consisting of boundary chords induced by the given normal
    surface surf.

    Consider the collection of boundary annuli given by cutting boundary tori
    along boundary curves of surf. The returned set will contain exactly one
    boundary chord spanning each such boundary annulus.

    The total number of returned boundary chords is therefore equal to the
    number of boundary components of surf. In particular, if surf is disjoint
    from the real boundary of its ambient triangulation, then this routine
    will return an empty set.

    For a boundary annulus containing a vertex of surf.triangulation(), the
    spanning boundary chord will be chosen to consist of two type-1 segments
    that meet at said vertex (in particular, the chord will have length 2).

    For any other boundary annulus, the spanning boundary chord is
    (necessarily) built entirely from a single type-2 segment.

    The spanning boundary chords will never have their ends abstractly joined
    to any other chords. Orientations will be chosen arbitrarily, subject to
    the constraint that boundary chords lying in the same boundary torus are
    oriented consistently with each other.

    Pre-condition:
    --> The given surf should be either a 2-sphere, projective plane, disc,
        annulus, or Mobius band. Moreover, if surf has a trivial boundary
        curve, then it must be a disc.
    --> The ambient triangulation surf.triangulation() should have minimal
        toroidal boundary.
    """
    if not surf.hasRealBoundary():
        return set()
    boundaryChords = set()
    for bc in surf.triangulation().boundaryComponents():
        # From the pre-conditions, bc is a two-triangle torus.
        bdryFace = bc.triangle(0)
        faceEmb = bdryFace.front()
        tet = faceEmb.tetrahedron()

        # Does the surface intersect bc? If so, then we might need to pick up
        # some boundary chords.
        normalArcs = [ surf.arcs( bdryFace.index(), v ).pythonValue()
                      for v in range(3) ]
        zeros = normalArcs.count(0)
        if zeros == 3:
            # Surface has no intersection with bc, so we can't pick up any
            # new boundary chords here.
            continue
        if zeros == 0:
            # Surface has a trivial boundary curve.
            #
            # From the pre-conditions, the surface must in fact be a disc,
            # which means that there cannot be any boundary chords at all.
            assert not boundaryChords
            return boundaryChords
        numBdryCurves = pythonGCD(*normalArcs)
        assert numBdryCurves <= 2, \
                "Failed pre-condition: Too many boundary curves"

        # Find boundary chord incident to the central faces in bc.
        if zeros == 2:
            # Let e denote the edge of bdryFace opposite the nonzero normal
            # arc. Take the spanning boundary chord to consist of two type-1
            # segments stradding one end (say, end 0) of e.
            v = normalArcs.index(numBdryCurves)
            oppEndpoints = {0,1,2,3} - { faceEmb.vertices()[v],
                                     faceEmb.vertices()[3] }
            oppEdge = tet.edge(*oppEndpoints)
            oppFront = oppEdge.front()
            oppBack = oppEdge.back()
            frontTet = oppFront.tetrahedron()
            backTet = oppBack.tetrahedron()
            frontEdgeEnds = [ oppFront.vertices()[0],
                             oppFront.vertices()[2] ]
            backEdgeEnds = [ oppBack.vertices()[0],
                            oppBack.vertices()[3] ]
            frontEdgeIndex = frontTet.edge(*frontEdgeEnds).index()
            backEdgeIndex = backTet.edge(*backEdgeEnds).index()

            # Arbitrarily choose to orient the chord from front to back
            # (i.e., orient front segment towards the vertex, and back
            # segment away from the vertex).
            frontEdgeMapping = frontTet.edgeMapping(
                    Edge3.faceNumber(*frontEdgeEnds) )
            backEdgeMapping = backTet.edgeMapping(
                    Edge3.faceNumber(*backEdgeEnds) )
            if frontEdgeMapping[0] == oppFront.vertices()[0]:
                frontSegPos = 0
                frontOrientation = -1
            else:
                frontSegPos = surf.edgeWeight(
                        frontEdgeIndex ).pythonValue()
                frontOrientation = 1
            if backEdgeMapping[0] == oppBack.vertices()[0]:
                backSegPos = 0
                backOrientation = 1
            else:
                backSegPos = surf.edgeWeight(
                        backEdgeIndex ).pythonValue()
                backOrientation = -1
            frontSeg = OrientedSegment(
                    surf, frontEdgeIndex, frontSegPos, frontOrientation )
            backSeg = OrientedSegment(
                    surf, backEdgeIndex, backSegPos, backOrientation )
            centralBdryChord = NormalChord( [ frontSeg, backSeg ] )
        elif zeros == 1:
            # Let v denote the vertex of bdryFace at which we have the zero
            # normal arc. Take the spanning boundary chord to consist of two
            # type-1 segments straddling v.
            v = normalArcs.index(0)
            segmentsInChord = []
            firstSegment = True
            for other in range(3):
                if other == v:
                    continue
                endpoints = { faceEmb.vertices()[v],
                             faceEmb.vertices()[other] }
                segEdgeMapping = tet.edgeMapping(
                        Edge3.faceNumber(*endpoints) )
                segEdgeIndex = tet.edge(*endpoints).index()
                segEdgeWeight = surf.edgeWeight(segEdgeIndex).pythonValue()
                if firstSegment:
                    # Orient towards the vertex. This ensures that the
                    # orientation matches the direction of traversal through
                    # the chord.
                    if segEdgeMapping[0] == faceEmb.vertices()[v]:
                        segPos = 0
                        segOrientation = -1
                    else:
                        segPos = segEdgeWeight
                        segOrientation = 1

                    # The next segment we encounter (obviously) won't be the
                    # first.
                    firstSegment = False
                else:
                    # Orient away from the vertex.
                    if segEdgeMapping[0] == faceEmb.vertices()[v]:
                        segPos = 0
                        segOrientation = 1
                    else:
                        segPos = segEdgeWeight
                        segOrientation = -1
                segmentsInChord.append( OrientedSegment(
                    surf, segEdgeIndex, segPos, segOrientation ) )
            centralBdryChord = NormalChord(segmentsInChord)
        else:   # Impossible.
            raise AssertionError()
        boundaryChords.add(centralBdryChord)

        # If surf has two boundary curves in bc, then we also have another
        # boundary chord spanning a boundary annulus built entirely from
        # parallel faces. Such a boundary chord necessarily consists entirely
        # of a single type-2 segment.
        if numBdryCurves == 2:
            # To make it easy to fulfil the promise that the chord
            # orientations are consistent, we simply choose a segment in the
            # same edge as one of the segments already in centralBdryChord.
            centralSeg = centralBdryChord[0]
            parallelBdryChord = NormalChord( [ OrientedSegment(
                surf, centralSeg.edgeIndex(), 1, centralSeg.orientation() ) ] )
            boundaryChords.add(parallelBdryChord)

            # We have already found two boundary chords, which (from the
            # pre-conditions) is the maximum possible.
            break

    # Done!
    return boundaryChords


def _boundsOrbitalCompressionDisc(chord):
    """
    Does the given normal chord bound an orbital compression disc?
    """
    if len(chord) != 2:
        return False
    for seg in chord:
        if seg.segmentType() != 1:
            return False

    # We are looking at a normal chord built precisely from two type-1
    # segments. From the definition, this chord bounds an orbital compression
    # disc if and only if the two segments belong to the same type-1 orbit.
    mySeg, yourSeg = chord
    return ( mySeg.translateAlongParallelCells(
        { yourSeg.reversed() } ) is not None )


#TODO Replace all usage of this with the new routines above.
def newIdealLoopEmbs( surf, oldLoops=[] ):
    """
    Returns surviving edge embeddings which describe the ideal loops after
    crushing the given normal surface surf.

    The given oldLoops list (which may be empty, and is empty by default)
    should be a list of pre-existing ideal loops, encoded as instances of
    IdealLoop. Each of these ideal loops must lie in the same triangulation
    as surf, and these ideal loops must all be mutually disjoint.

    If there are no pre-existing ideal loops, then surf should be of one of
    the following types:
    --> A 2-sphere.
    --> A disc.
    --> An annulus with nontrivial boundary curves.
    --> A Mobius band with nontrivial boundary curve.
    Otherwise, letting W denote the weight of surf on the pre-existing ideal
    loops, surf should be of one of the following types:
    --> A 2-sphere with either W == 2 or W == 0.
    --> A disc with W == 1 and nontrivial boundary curve.
    --> A disc with W == 0.
    --> A projective plane with W == 1.
    This routine raises ValueError if surf is not of one of these allowed
    types.

    We also require surf to be a quadrilateral vertex normal surface, but
    this routine does not check this condition.

    This routine returns a list describing the ideal loops that would arise
    after crushing the given surface (see below for a more detailed
    description of how the ideal loops before crushing are related to the
    ideal loops after crushing). Each such ideal loop is encoded as a list of
    surviving edge embeddings.

    A caveat to this is that when the given surf is a 2-sphere, there is one
    possible degenerate ideal loop: a pair of edges giving an unknotted loop,
    such that the two edges get merged to become a single non-loop edge after
    crushing. This routine does not check for such degenerate loops, so they
    might appear in the returned list.

    Crushing the given surface has the following effects:
    --> Pre-existing ideal loops that are disjoint from the surface will be
        left topologically untouched. In particular, their orientations will
        be preserved.
    --> Ideal loops that intersect the surface will be split into multiple
        chords, and each such chord may or may not survive to become a new
        ideal loop after crushing. The orientation will be preserved for the
        chords that do survive.
    --> If the surface is an annulus (which, as specified above, must be
        disjoint from all pre-existing ideal loops), then crushing might
        create an entirely new ideal loop. This new loop will be assigned an
        arbitrary orientation.

    Pre-condition:
    --> The given surf should be a quadrilateral vertex normal surface.
    --> If surf is an annulus, then each boundary component that it meets
        must be a two-triangle torus.
    """
    # The given surf must be either a 2-sphere or an annulus. Moreover:
    # - In the 2-sphere case, we allow one of the ideal loops to have
    #   nonempty intersection with the surface.
    # - In the annulus case, we might create a new ideal loop by flattening
    #   a chain of boundary bigon faces.
    if isSphere(surf):
        loopMustBeDisjoint = False
        possibleLoopFromBoundary = False
    elif isAnnulus(surf):
        loopMustBeDisjoint = True
        possibleLoopFromBoundary = True
    else:
        allowed = "annuli and 2-spheres"
        msg = ( "This routine currently only accepts {} ".format(allowed) +
                "for the input surface." )
        raise ValueError(msg)

    # Find the ideal loops that arise from the pre-existing ideal loops.
    tri = surf.triangulation()
    newLoopEmbs = []
    survivors = OrientedSegment.survivors(surf)
    for oldLoop in oldLoops:
        wt = oldLoop.weight(surf)
        if wt == 2:
            if loopMustBeDisjoint:
                msg = ( "Too many intersections between the surface and " +
                        "the pre-existing ideal loops." )
                raise ValueError(msg)

            # For a 2-sphere, we currently only allow at most one ideal loop
            # to intersect the surface.
            loopMustBeDisjoint = True
        elif wt != 0:
            msg = ( "Each ideal loop must intersect the surface in " +
                    "either exactly 0 points or exactly 2 points." )
            raise ValueError(msg)

        # The given surface splits the current oldLoop into some number of
        # chords. Which of these chords survive to become new ideal loops
        # after crushing?
        for chord in oldLoop.splitIntoChords(surf):
            seg = chord[0]
            survivingSeg = seg.translateAlongParallelCells(survivors)
            if survivingSeg is None:
                # This chord does not survive after crushing.
                continue

            # This chord survives after crushing.
            newLoop = [ survivingSeg.survivingEmbedding() ]
            for seg in chord[1:]:
                survivingSeg = seg.translateAlongParallelCells(survivors)
                newLoop.append( survivingSeg.survivingEmbedding() )
            newLoopEmbs.append(newLoop)

    # Will there also be an entirely new ideal loop created by flattening a
    # chain of boundary bigons?
    if possibleLoopFromBoundary and countIncidentBoundaries(surf) == 1:
        # Find a segment incident to the chain of boundary bigons.
        for e in tri.edges():
            ei = e.index()
            if ( e.isBoundary() and
                    surf.edgeWeight(ei).pythonValue() >= 2 ):
                # Arbitrarily assign orientation +1.
                seg = OrientedSegment( surf, ei, 1, 1 )
                break

        # If this segment survives after crushing, then it will form a new
        # ideal loop of length one.
        survivingSeg = seg.translateAlongParallelCells(survivors)
        if survivingSeg is not None:
            newLoopEmbs.append( [ survivingSeg.survivingEmbedding() ] )

    # Done!
    return newLoopEmbs


#TODO Replace all usage of this with the new routines above.
def fillIdealEdges( tri, endpoints ):
    """
    Fills in two-triangle boundary 2-spheres of tri that have a symmetric pair
    of vertices in the given set of endpoints, and returns a list of the
    resulting ideal edges.

    The endpoints set should be specified as a set of integers corresponding
    to vertex indices in the given triangulation tri.

    In detail, each two-triangle boundary 2-sphere B that we fill in falls
    under one of the following two cases:
    --> If B is isomorphic to the boundary of a triangular pillow, then any
        two out of the three vertices of B form a symmetric pair. Let e denote
        the edge of B that has both vertices in the given endpoints set. We
        fill in B by gluing the two faces of B together in the obvious way.
        The edge e becomes one of the ideal edges in the returned list.
    --> If B is isomorphic to the boundary of a snapped 3-ball, then the only
        symmetric pair of vertices is given by the two "poles" (that is, the
        two vertices not incident to the equator edge). We fill in B by
        attaching a snapped 3-ball. This introduces a new edge connecting the
        two poles, which becomes one of the ideal edges in the returned list.

    This routine modifies the given triangulation directly. If the
    triangulation is currently oriented, then the filling operation will
    preserve this orientation.
    """
    fillEdges = []  # Items will be triples ( tet, edgeNum, doLayer ).
    for bc in tri.boundaryComponents():
        # Is bc a real two-triangle 2-sphere boundary component?
        #
        # Note that bc might have vertices pinched together (tri is not
        # assumed to be valid), so we must build the Triangulation2 to be sure
        # of computing the Euler characteristic correctly.
        built = bc.build()
        if built.eulerChar() != 2:
            continue
        if not bc.isReal():
            continue
        if bc.size() != 2:
            continue

        # Do we have a pillow-boundary or a snapped-ball-boundary?
        equatorEdgeNum = None
        for e in range(3):
            if built.triangle(0).adjacentTriangle(e) == built.triangle(1):
                if equatorEdgeNum is None:
                    equatorEdgeNum = e
                else:
                    # For snapped-ball-boundary, the equator edge is the
                    # unique edge that the two boundary triangles have in
                    # common. Since we have now found two common edges, we
                    # cannot have snapped-ball-boundary.
                    equatorEdgeNum = None
                    break

        # Do we have a symmetric pair of vertices in the endpoints set?
        if equatorEdgeNum is None:
            # We have pillow-boundary, so the answer is "yes" if and only if
            # exactly two out of the three vertices are in the endpoints set.
            idealEdgeNum = None
            for e in range(3):
                # In bc.triangle(0), the ideal edge should be opposite the
                # unique vertex that is *not* in the endpoints set.
                if bc.triangle(0).vertex(e).index() not in endpoints:
                    if idealEdgeNum is None:
                        idealEdgeNum = e
                    else:
                        idealEdgeNum = None
                        break
            if idealEdgeNum is not None:
                # We already have the ideal edge, so doLayer should be False.
                emb = bc.triangle(0).edge(e).front()
                fillEdges.append(
                        ( emb.tetrahedron(), emb.face(), False ) )
        else:
            # We have snapped-ball-boundary, so the answer is "yes" if and
            # only if the "poles" (ie, the vertices not incident to the
            # equator edge) both lie in the endpoints set, and the vertex on
            # the equator edge does not lie in the endpoints set.
            if ( bc.triangle(0).edge(
                equatorEdgeNum ).vertex(0).index() in endpoints ):
                # Equator vertex is in the endpoints set.
                # Move on to the next boundary component.
                continue
            if bc.triangle(0).vertex(equatorEdgeNum).index() not in endpoints:
                # Pole vertex is not in the endpoints set.
                # Move on to the next boundary component.
                continue

            # We just need to check the pole incident to bc.triangle(1) now.
            otherPole = built.triangle(0).adjacentEdge(equatorEdgeNum)
            if bc.triangle(1).vertex(otherPole).index() in endpoints:
                # Later, we will need to layer across the equator edge to
                # obtain the ideal edge (ie, we need doLayer to be True).
                emb = bc.triangle(0).edge(equatorEdgeNum).front()
                fillEdges.append(
                        ( emb.tetrahedron(), emb.edge(), True ) )

    # Perform all the fillings.
    idealEdgeLocations = []
    for tet, edgeNum, doLayer in fillEdges:
        if doLayer:
            # Make sure to use oriented layering.
            tet = layerOn( tet.edge(edgeNum) )
            edgeNum = 5
        idealEdgeLocations.append( ( tet, edgeNum ) )
        idealEdge = tet.edge(edgeNum)

        # Close up the pillow boundary.
        #
        #       back    front
        #          0    0
        #          •    •
        #         /|    |\
        #        / |    | \
        #      3•  |    |  •2
        #        \ |    | /
        #         \|    |/
        #          •    •
        #          1    1
        #
        front = idealEdge.front()
        back = idealEdge.back()
        front.tetrahedron().join(
                front.vertices()[3],
                back.tetrahedron(),
                back.vertices() * Perm4(2,3) * front.vertices().inverse() )

    # All done!
    return [ tet.edge(edgeNum) for tet, edgeNum in idealEdgeLocations ]
