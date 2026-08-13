"""
Find the ideal edges after crushing a normal surface.
"""
from enum import Enum, auto
from math import gcd as pythonGCD
from regina import *
from loop import IdealLoop  #TODO Probably need to keep this, but double-check.
from triloops import EdgeIdealTriangulation
from chord import pairUpChordEndsByCrushing
from segment import OrientedSegment
from aux.bdry import hasOnlyMinimalRealTorusBoundaryComponents
from aux.tetrenum import tetRenumbering
from aux.quad import tetHasQuads
from aux.surface import SurfaceType, hasOnlyNonTrivialBoundaryCurves
from aux.surface import isSphere, isAnnulus, countIncidentBoundaries
from aux.looperror import NotLoop
from retriangulate.insert import layerOn


def checkCrushAllowedInRealTriangulation(surf):
    """
    Checks that we are allowed to crush along the given normal surface.

    This routine raises ValueError if and only if
    realTriangulationAllowsCrush(surf) returns False.

    Pre-condition:
    --> surf.triangulation() is orientable().
    """
    if surfType == SurfaceType.ANNULUS:
        if hasOnlyNonTrivialBoundaryCurves(surf):
            return
        raise ValueError( "To crush along an annulus, both of its " +
                         "boundary curves must be nontrivial" )
    elif surfType == SurfaceType.MOBIUS:
        if hasOnlyNonTrivialBoundaryCurves(surf):
            return
        raise ValueError( "To crush along a Mobius band, its boundary " +
                         "curve must be nontrivial" )
    elif surfType in { SurfaceType.SPHERE, SurfaceType.DISC }:
        return
    raise ValueError( "With a real triangulation, we cannot crush along a " +
                     "surface of type {}".format(surfType.name) )


def realTriangulationAllowsCrush(surf):
    """
    Are we allowed to crush along the given normal surface?

    This routine returns True if and only if surf is of one of the following
    types:
    --> A 2-sphere.
    --> A disc.
    --> An annulus with nontrivial boundary curves.
    --> A Mobius band with nontrivial boundary curve.

    Pre-condition:
    --> surf.triangulation() is orientable().
    """
    try:
        checkCrushAllowedInRealTriangulation(surf)
    except ValueError:
        return False
    return True


#TODO Overhaul decomposeAlong() to not only return the components that
#       survive crushing, but also to do the book-keeping of tracking deleted
#       components and counting the number of orbital compressions.
#TODO Update usage for new output format.
def decomposeAlong( surf, edgeIdealTri=None ):
    """
    Decomposes along surf, and returns a list of the resulting components.

    In detail, each item in the returned list will be instance of either
    EdgeIdealTriangulation or Triangulation3. The underlying triangulations
    will just be the components that result from crushing surf.triangulation()
    along surf. If there are any EdgeIdealTriangulation objects in the
    returned list, then the IdealLoop objects that are embedded in these
    triangulations indicate loops that should be drilled out to obtain
    topologically useful 3-manifolds obtained by decomposing along surf.

    If there are no pre-existing IdealLoop objects to track, then only surf
    needs to be supplied as input to this routine. Otherwise, both surf and
    edgeIdealTri should be supplied, in which case edgeIdealTri should be an
    instance of EdgeIdealTriangulation that tracks all of the pre-existing
    IdealLoop objects.

    This routine might (but is not guaranteed to) detect that some ideal loops
    are "trivial" in the sense that they bound embedded discs. This routine
    handles such a trivial ideal loop L in one of two ways:
    --> If L is the only ideal loop in its ambient triangulation, then this
        routine will silently delete L.
    --> Otherwise, if there are other ideal loops in the ambient triangulation
        of L, then this routine will raise BoundsDisc.

    A side-effect of this routine is that it might (but is not guaranteed to)
    detect and silently delete some (or all, or none) of the ideal loops that
    are "trivial" in the sense that they bound embedded discs.

    If edgeIdealTri is None, then realTriangulationAllowsCrush(surf) must be
    True. Otherwise, edgeIdealTri.allowsCrush(surf) must be True. This
    routine raises ValueError if surf does not satisfy these conditions.

    We also require surf to be a quadrilateral vertex normal surface, but
    this routine does not check this condition.

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

    #TODO Make new buildNewLoopsFromIdealChords() and
    #       buildNewLoopsFromBoundaryChords() routines.
    #
    #       There's still a question about whether extracting surviving edge
    #       embeddings should happen inside or outside these new routines.
    #
    #       Possibly the answer is to do this *outside*, provided we can
    #       detect orbital compressions *inside*.

    # Check that surf is of one of the required types.
    surfType = SurfaceType.recognise(surf)
    if edgeIdealTri is None:
        checkCrushAllowedInRealTriangulation(surf)
    else:
        # Enforce the precondition that the two input objects reference
        # precisely the same Triangulation3 object in memory.
        if tri is not edgeIdealTri.triangulation():
            raise RuntimeError( "decomposeAlong() requires the input " +
                               "NormalSurface and the input " +
                               "EdgeIdealTriangulation to reference the " +
                               "same Triangulation3 object in memory" )
        edgeIdealTri.checkCrushAllowed(surf)

    #TODO Decide whether to bother with allowing projective planes as input.

    #TODO Make a final decision on how to deal with trivial loops, and then
    #       check that trivial loops are actually dealt with as intended.

    #TODO Don't forget that whenever we close up an invalid boundary
    #       2-sphere, our length-2 boundary chord construction ensures that
    #       we can always shorten/redirect/whatever our ideal loop along the
    #       newly closed-up triangular face.

    # Find where the new ideal loops will be after crushing.
    #
    #NOTE The newIdealLoopEmbs() requires, and will check, that the given
    #       surface is of one of the allowed types.
    loopEmbSeqsInOldTri = newIdealLoopEmbs( surf, oldLoops )
    doomed = [ tet for tet in surf.triangulation().tetrahedra()
              if tetHasQuads( surf, tet.index() ) ]
    tetIndicesAfterCrush = tetRenumbering(doomed)
    crushed = surf.crush()
    loopEmbSeqs = []
    for oldEmbSequence in loopEmbSeqsInOldTri:
        embSequence = []
        for oldEmb in oldEmbSequence:
            crushedTet = crushed.tetrahedron(
                    tetIndicesAfterCrush[ oldEmb.tetrahedron().index() ] )
            embSequence.append( EdgeEmbedding3(
                crushedTet, oldEmb.vertices() ) )
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
    output = []
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

            # Note that we could have a degenerate loop.
            try:
                loop = IdealLoop( edgeList, orientation )
            except NotLoop:
                #TODO Is ignoring these the right thing to do?

                # Ignore degenerate loop.
                continue
            else:
                loops.append(loop)

        # If we have any loops at all, then package them all together as a
        # single EdgeIdealTriangulation. Otherwise, just add an ordinary
        # Triangulation3 to the output list.
        #
        #TODO Consider simplifying before adding to the output list. But this
        #       requires deciding how this routine should behave if
        #       simplification raises BoundsDisc.
        if loops:
            output.append( EdgeIdealTriangulation(loops) )
        else:
            output.append(tri)

    # All done!
    return output


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


#TODO Still need to make a final decision on the name for this routine.
def findNewIdealLoops( surf, edgeIdealTri ):
    """
    Finds the new ideal loops that would arise from crushing the given normal
    surface surf.

    This routine returns a list describing the new ideal loops. In detail:
    --> Pre-existing ideal loops in edgeIdealTri that are disjoint from the
        surface will be left topologically untouched. In particular, their
        orientations will be preserved.
    --> Ideal loops that intersect the surface will be split into multiple
        ideal chords; additionally, if the surface is incident to the
        boundary, there might be some boundary chords. Each such (ideal or
        boundary) chord either survives the crushing operation, or is
        entirely destroyed by crushing. For the surviving chords, crushing
        will essentially rearrange how the endpoints of these chords are
        joined together, thereby yielding new ideal loops.
    Each element of the returned list describes one of the new ideal loops
    via a pair consisting of the following items:
    (0) A (possibly empty) list of surviving edge embeddings, appearing in
        order of traversal around the ideal loop, and also oriented
        consistently with the order of traversal.
    (1) A status given by _IdealLoopStatus.

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
    newLoops = []   # We populate and return this list.
    survivors = OrientedSegment.survivors(surf)
    #TODO Do we allow edgeIdealTri to be None in this function, or handle
    #       that case in a separate function?
    internalChords = edgeIdealTri.splitIntoChords(surf)
    boundaryChords = _findBoundaryChords(surf)

    # Extract all the segments from the internal chords and put them together
    # to form the new loops. If one of the internal chords has unjoined ends,
    # then these will need to be joined to the boundary chord (which is
    # unique if it exists, since the pre-conditions ensure that the only case
    # where we have a boundary chord is when surf is a disc).
    #
    # This also checks for inconsistently oriented chords.
    chordSequences = []
    while internalChords:
        currentChord = internalChords.pop()
        chordsInNewLoop = [currentChord]
        loopStatus = _IdealLoopStatus.CONSISTENT    # Until proven otherwise.
        currentTailEnd = currentChord.joinedEnd(1)
        currentChord = currentChord.joinedChord(1)
        if currentChord is None:
            # From the pre-conditions, surf must be a disc with edge-ideal
            # weight 1, and we must have found the unique internal chord
            # which needs to be joined with the unique boundary chord.
            assert len(boundaryChords) == 1
            bdryChord = boundaryChords.pop()
            pairUpChordEndsByCrushing( chordsInNewLoop[0], bdryChord )

            # The current choice of orientation on the boundary chord is
            # arbitrary and meaningless. We re-orient (if necessary) to
            # ensure that the orientation is consistent with the pre-existing
            # internal chord.
            if chordsInNewLoop[0].joinedEnd(0) == 0:
                chordsInNewLoop.append( bdryChord.reversed() )
            else:
                chordsInNewLoop.append(bdryChord)
        else:
            # Traverse the new loop, and pick up all its constituent internal
            # chords.
            while currentChord != chordsInNewLoop[0]:
                internalChords.remove(currentChord)
                if currentTailEnd == 1:
                    # The current chord is oriented inconsistently with the
                    # first chord in this new loop.
                    loopStatus = _IdealLoopStatus.INCONSISTENT
                    chordsInNewLoop.append( currentChord.reversed() )
                else:
                    chordsInNewLoop.append(currentChord)
                currentHeadEnd = 1 - currentTailEnd

                # Move on to the next chord in the loop.
                currentTailEnd = currentChord.joinedEnd(currentTailEnd)
                currentChord = currentChord.joinedChord(currentTailEnd)
            assert currentTailEnd == 0  # Tail of the first chord
        chordSequences.append( ( chordsInNewLoop, loopStatus ) )

    #TODO What about chords that still remain in the boundaryChords set? If
    #       edgeIdealTri is not None, then I think we can just throw these
    #       away (this is the case where we crushed a disc with edge-ideal
    #       weight 0, and hence don't obtain a new ideal loop from the
    #       (unique) boundary chord).

    # Extract surviving embeddings (if any).
    #
    # This also checks for orbital compressions.
    newLoops = []
    for chordsInNewLoop, loopStatus in chordSequences:
        # Check whether this loop is compressed away by an orbital
        # compression disc.
        #
        # For detecting orbital compression discs bounded by real (as opposed
        # to edge-ideal) boundary edges, this check is sufficient because of
        # the promise made in _findBoundaryChords() that we use chords of
        # length 2 whenever possible.
        if len(chordsInNewLoop) == 1 and len( chordsInNewLoop[0] ) == 2:
            assert loopStatus == _IdealLoopStatus.CONSISTENT
            mySeg, yourSeg = *chordsInNewLoop[0]
            if ( mySeg.translateAlongParallelCells(
                { yourSeg.reversed() } ) is not None ):
                # We have an orbital compression disc here.
                loopStatus = _IdealLoopStatus.COMPRESSED

        # As a consequence of the pre-condition that surf is quad vertex,
        # every orbit is either simple, or deformation retracts to an orbital
        # compression disc. With not much work, it follows that one segment
        # of the current loop survives if and only if every segment in every
        # chord of the current loop survives.
        newLoopEmbeddings = []
        newLoopSurvives = True  # Until proven otherwise.
        for chord in chordsInNewLoop:
            for seg in chord:
                survivingSeg = seg.translateAlongParallelCells(survivors)
                if survivingSeg is None:
                    # The current new loop doesn't survive, so we immediately
                    # break out and add an empty collection of surviving
                    # embeddings to the newLoops list.
                    newLoopSurvives = False
                    break
                else:
                    newLoopEmbeddings.append(
                            survivingSeg.survivingEmbedding() )
            if not newLoopSurvives:
                break
        newLoops.append( ( newLoopEmbeddings, loopStatus ) )

    # All done!
    return newLoops


def _findBoundaryChords(surf):
    """
    Returns a set consisting of boundary chords induced by the given normal
    surface.

    Consider the collection of boundary annuli given by cutting boundary tori
    along boundary curves of surf. The returned set will contain exactly one
    boundary chord spanning each such boundary annulus. For boundary annuli
    containing a vertex of surf.triangulation(), it is guaranteed that the
    corresponding boundary chord will be chosen to consist of two type-1
    segments (in other words, the chord will have length 2).

    Observe that the total number of returned boundary chords is therefore
    equal to the number of boundary components of surf.

    The ends of the returned chords are never abstractly joined to any other
    chords. The orientations on the returned chords are chosen arbitrarily.

    In the particular case where surf is disjoint from the real boundary of
    its ambient triangulation, this routine returns the empty set.

    Pre-condition:
    --> The given surf should be either a 2-sphere, projective plane, disc,
        annulus, or Mobius band. Moreover, if surf has a trivial boundary
        curve, then it must be a disc.
    --> The ambient triangulation surf.triangulation() should have minimal
        toroidal boundary.
    """
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
        numBdryCurves = pythonGCD(*normalArcs)
        assert numBdryCurves <= 2, \
                "Failed pre-condition: Too many boundary curves"
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

        # Find possible boundary chord sandwiched between two parallel
        # boundary curves.
        if numBdryCurves == 2:
            # A boundary edge opposite a zero normal arc coordinate will
            # always have positive weight.
            v = normalArcs.index(0)
            endpoints = {0,1,2,3} - { faceEmb.vertices()[v],
                                     faceEmb.vertices()[3] }
            oppEdgeIndex = tet.edge(*endpoints).index()
            segPos = 1  # Any odd segment position will do.

            # Arbitrarily choose the orientation to be +1.
            boundaryChords.add( NormalChord( [ OrientedSegment(
                surf, oppEdgeIndex, segPos, 1 ) ] ) )

        # Find possible boundary chord incident to the central faces in bc.
        if zeros == 2:
            # Take the boundary chord to consist of two type-1 segments
            # straddling one end (say, end 0) of the edge of bdryFace
            # opposite the *nonzero* normal arc.
            v = normalArcs.index(numBdryCurves)
            endpoints = {0,1,2,3} - { faceEmb.vertices()[v],
                                     faceEmb.vertices()[3] }
            oppEdge = tet.edge(*endpoints)
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
            frontEdgeMapping = frontTet.edgeMapping(*frontEdgeEnds)
            backEdgeMapping = backTet.edgeMapping(*backEdgeEnds)
            if frontEdgeMapping[0] == oppFront.vertices[0]:
                frontSegPos = 0
                frontOrientation = -1
            else:
                frontSegPos = surf.edgeWeight(
                        frontEdgeIndex ).pythonValue()
                frontOrientation = 1
            if backEdgeMapping[0] == oppBack.vertices[0]:
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
            boundaryChords.add( NormalChord( [ frontSeg, backSeg ] ) )
        elif zeros == 1:
            # Take the boundary chord to consist of two type-1 segments
            # straddling the vertex at which we have the *zero* normal arc.
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
            boundaryChords.add( NormalChord(segmentsInChord) )
        else:   # Impossible.
            raise AssertionError()

    # Done!
    return boundaryChords


#TODO How best to deal with ideal chords which are joined to each other, and
#       with ideal chords arising from closing up invalid vertices?
def idealLoopsFromRealBoundary(surf):
    """

    This routine returns a list describing the new ideal chords that would
    arise from the real boundary after crushing the given normal surface
    surf.
    ... TODO ...
    Each element of the returned list describes one of the new ideal chords
    via a surviving edge embedding.
    """
    #TODO Find new loops arising from real boundary.


#TODO Update documentation and implementation to:
#       --> use the new TriangulationWithEmbeddedLoops class, and
#       --> account for the extra cases that arise from SFS recognition.
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
    #TODO Update the checks below to allow for the extra cases that arise in
    #       the context of bounded SFS recognition.

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

    #TODO If we crushed an annulus, it would probably be useful to use
    #   fillIdealEdges() to include additional ideal loops obtained by filling
    #   in pinched 2-sphere boundary components.
    #
    #   If/when we implement this functionality, we will need to document the
    #   possibility that we could create an additional new ideal loop. We
    #   should probably also note that this would come at the cost of
    #   introducing a new tetrahedron.

    # Done!
    return newLoopEmbs


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
