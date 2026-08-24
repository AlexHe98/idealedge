"""
Scratch work for ideal edges.
"""
from sys import argv
from timeit import default_timer
from regina import *
from idealedge import decomposeAlong, newIdealLoopEmbs, fillIdealEdges
from idealedge import ComponentDeletedByCrushing as DelComp
from idealedge import SurfaceToCrushInSuspectedSFS as CandidateSurface
from loop import IdealLoop, BoundsDisc
from triloops import EdgeIdealTriangulation
from pinch import drillMeridian
from wedge import NonSurvivingTriangularOrbitType as OrbitType
from wedge import nonSurvivingTriangularOrbitCounts as orbitCounts
from construct.sfs import orientableSFS
from aux.tetrenum import tetRenumbering
from aux.quad import tetHasQuads
from aux.surface import isSphere, isAnnulus


#TODO Make this more general.
def filledHomology(annulus):
    """
    Given an annulus A in a triangulation whose boundary is a two-triangle
    torus B, such that both boundary curves of A are parallel to an edge of B,
    computes the first homology of the manifold given by filling along the
    slope of the boundary of A.
    """
    # What's the homology of the manifold before filling?
    tri = annulus.triangulation()
    markedH1 = tri.markedHomology(1)
    rank = markedH1.rank()
    numInvFacs = markedH1.countInvariantFactors()
    invFacs = [ markedH1.invariantFactor(i) for i in range(numInvFacs) ]
    numGenerators = rank + numInvFacs
    # Each column represents a generator.
    # Each row gives a relation.
    presentation = []
    for i in range(numInvFacs):
        relation = [0] * numGenerators
        relation[i] = invFacs[i]
        presentation.append(relation)

    # We need to add one extra relation corresponding to the curve c that is
    # killed by the Dehn filling. By assumption, there is a boundary edge that
    # is parallel to c.
    missedEdgeIndex = None
    for e in tri.edges():
        if not e.isBoundary():
            continue
        if annulus.edgeWeight( e.index() ).pythonValue() == 0:
            missedEdgeIndex = e.index()
            break
    if missedEdgeIndex is None:
        raise ValueError( "Given surface doesn't satisfy preconditions" )

    # Write the missed edge in terms of generators, so that we can then add a
    # relation that kills c.
    cycle = [0] * tri.countEdges()
    cycle[missedEdgeIndex] = 1
    relation = markedH1.snfRep(cycle)
    #TODO BEGIN TEST
    if presentation:
        print( AbelianGroup( MatrixInt(presentation) ) )
    print(relation)
    #TODO END TEST
    presentation.append(relation)
    return AbelianGroup( MatrixInt(presentation) )


def crushCandidateVerticalSurfaces( surfaces, threshold=30 ):
    """
    Crushes all candidate vertical surfaces in the given list of quad vertex
    normal surfaces.

    This routine prints details of the crushed triangulations.

    This routine attempts to identify the topology of the manifold that
    results from crushing. The main strategy is to simplify and attempt
    combinatorial recognition. Additionally, whenever this routine encounters
    a component whose number of tetrahedra is strictly less than the
    threshold (default 30), it will also use more computationally intensive
    recognition algorithms involving normal surfaces.

    Pre-condition:
    --> Every boundary component of the ambient triangulation must be a
        two-triangle torus.
    """
    start = default_timer()
    annulusCount = 0
    for surfNum, surf in enumerate(surfaces):
        if CandidateSurface.recognise(surf) != CandidateSurface.VERTICAL:
            continue
        thin = surf.isThinEdgeLink()
        if thin[0] is not None:
            # Don't bother with thin edge links.
            continue
        annulusCount += 1
        print()
        print( "Time: {:.6f}. Crush #{}.".format(
            default_timer() - start, surfNum) )

        # Crush this candidate vertical surface.
        decomposed, numOrbCuts, delComps, inconsistent = decomposeAlong(surf)
        twists = []
        for _ in range( delComps[ DelComp.FIBRE_PLUS ] ):
            twists.append(1)
        for _ in range( delComps[ DelComp.FIBRE_MINUS ] ):
            twists.append(-1)
        crushedDesc = "Crushed #{}".format(surfNum)
        if not surf.isOrientable():
            crushedDesc += " Lost (2,1)."
        for twist in twists:
            crushedDesc += " Lost (3,{}).".format(twist)
        print(crushedDesc)
        if inconsistent:
            print( "--------------------" )
            print( "NON-ORIENTABLE BASE!" )
            print( "--------------------" )
        #TODO Use numOrbCuts (number of orbital compressions).
        #TODO Update filledHomology() and use it for sanity checking.
        for edgeIdealTri in decomposed:
            if isinstance( edgeIdealTri, EdgeIdealTriangulation ):
                try:
                    #TODO TEST
                    edgeIdealTri.simplifyMonotonic()
#                    edgeIdealTri.simplify()
                except BoundsDisc:
                    #TODO
                    print( "Loop bounds disc!" )
                else:
                    # Try to identify the topology of edgeIdealTri.
                    decomposedDesc = "Decomposed into fibres: "
                    decomposedList = decomposeAlongSpheres(edgeIdealTri)
                    for newEdgeIdealTri in decomposedList:
                        if len(newEdgeIdealTri) != 1:
                            decomposedDesc +=\
                                    "N/A (piece w/ {} loops)".format(
                                            len(newEdgeIdealTri) )
                            continue

                        # We should be able to drill and get Seifert fibre
                        # parameters.
                        try:
                            newTriWithMeridian = drillMeridian(
                                    newEdgeIdealTri[0] )
                        except BoundsDisc:
                            decomposedDesc += "N/A (S2 x S1); "
                        else:
                            newTriWithMeridian.minimiseBoundary()
                            newTriWithMeridian.simplify()
                            newTriWithMeridian.simplify()

                            # There is only one BoundaryLoop, corresponding
                            # to the meridian. Also, because we minimised the
                            # boundary, the meridian is guaranteed to be
                            # given by a single edge.
                            newMerEdgeIndex = newTriWithMeridian[0][0]
                            newSurf = newTriWithMeridian.triangulation().nonTrivialSphereOrDisc()
                            if newSurf is None:
                                decomposedDesc += "N/A (no disc); "
                            elif newSurf.eulerChar() == 2:
                                decomposedDesc += "unknown (found sphere); "
                            else:
                                decomposedDesc += "({},{}); ".format(
                                        *fibreParams( newSurf, newMerEdgeIndex ) )

                    # Format decomposedDesc correctly.
                    decomposedDesc = decomposedDesc[:-2]
                    print( "    " + decomposedDesc )
                # End of try-except-else
            else:
                print( "Component with no loops!" )
        # End of loop through decomposed list.

    # All done!
    print()
    print( "Time: {:.6f}. All done!".format(
        default_timer() - start ) )
    return


def fibreParams( surf, merEdgeIndex ):
    # Use boundary edge weights of the disc to calculate
    # Seifert parameters.
    drilled = surf.triangulation()
    merWt = surf.edgeWeight(merEdgeIndex).pythonValue()
    merEdge = drilled.edge(merEdgeIndex)
    front = merEdge.front()
    ver = front.vertices()
    tet = front.tetrahedron()
    lower = tet.edge( ver[0], ver[2] )
    upper = tet.edge( ver[1], ver[2] )
    lowWt = surf.edgeWeight( lower.index() ).pythonValue()
    uppWt = surf.edgeWeight( upper.index() ).pythonValue()
    if merWt == lowWt + uppWt:
        print("M=L+U")
        shift = lowWt
    elif uppWt == merWt + lowWt:
        print("U=M+L")
        shift = -lowWt
    elif lowWt == merWt + uppWt:
        print("L=M+U")
        shift = uppWt
    else:
        raise ValueError( "Weights don't add up." )
    #q = pow( shift, -1, merWt )
    q = shift % merWt
    if q > merWt // 2:
        q -= merWt
    return ( merWt, q )


def decomposeAlongSpheres(edgeIdealTri):
    """
    Returns a list of ideal loops given by repeatedly decomposing the given
    ideal loop along spheres that intersect the ideal loop twice.
    """
    # Repeatedly crush along spheres that intersect the ideal loop at most
    # twice.
    toProcess = [edgeIdealTri]
    decomposedList = []
    while toProcess:
        oldEdgeIdealTri = toProcess.pop()
        tri = oldEdgeIdealTri.triangulation()

        # Search for a suitable sphere to crush.
        enumeration = TreeEnumeration( tri, NormalCoords.Quad )
        while True:
            # Get the next 2-sphere.
            if enumeration.next():
                sphere = enumeration.buildSurface()
                if not isSphere(sphere):
                    continue
            else:
                # No suitable 2-sphere means that we're done with the current
                # oldEdgeIdealTri.
                decomposedList.append(oldEdgeIdealTri)
                break

            # Does the sphere intersect the ideal loop at most twice?
            wt = oldEdgeIdealTri.weight(sphere)
            if wt != 2:
                #TODO Actually do something with the following cases.
                if wt == 0:
                    print( "Found sphere disjoint from ideal loop!" )
                if wt == 1:
                    print( "Found sphere intersecting ideal loop once!" )

                # Continue searching for suitable spheres.
                continue

            # See what happens if we crush.
            decomposed, numOrbCuts, delComps, inconsistent = decomposeAlong(
                    sphere, oldEdgeIdealTri )
            twists = []
            for _ in range( delComps[ DelComp.FIBRE_PLUS ] ):
                twists.append(1)
            for _ in range( delComps[ DelComp.FIBRE_MINUS ] ):
                twists.append(-1)
            lostFibres = ""
            for twist in twists:
                lostFibres += " Lost (3,{}).".format(twist)
            if lostFibres:
                print( lostFibres[1:] )
            if inconsistent:
                print( "vvvvvvvvvvvvvvvvvvvv" )
                print( "NON-ORIENTABLE BASE!" )
                print( "^^^^^^^^^^^^^^^^^^^^" )
            #TODO Use numOrbCuts (number of orbital compressions).
            for newEdgeIdealTri in decomposed:
                if isinstance( newEdgeIdealTri, EdgeIdealTriangulation ):
                    try:
                        #TODO TEST
                        newEdgeIdealTri.simplifyMonotonic()
#                        newEdgeIdealTri.simplify()
                    except BoundsDisc:
                        #TODO
                        print( "Loop bounds disc!" )
                    else:
                        toProcess.append(newEdgeIdealTri)

            # Found and crushed a suitable sphere, so stop enumerating.
            break

    # If we reach this point, then we have decomposed as far as possible, and
    # everything remaining has no suitable spheres.
    return decomposedList


if __name__ == "__main__":
    genus = int( argv[1] )
    boundaries = int( argv[2] )
    params = [ int(n) for n in argv[3:] ]
    fibres = []
    while params:
        q = params.pop()
        p = params.pop()
        fibres.append( (p,q) )
    print( "g={}, b={}".format( genus, boundaries ) )
    print(fibres)
    print()
    tri = orientableSFS( genus, boundaries, *fibres )
    tri.simplify()
    tri.simplify()
    #NOTE As of Regina 7.4, NS_QUAD and NS_VERTEX have been deprecated and
    #       replaced with NormalCoords.Quad and NormalList.Vertex,
    #       respectively.
    surfaces = NormalSurfaces( tri, NormalCoords.Quad, NormalList.Vertex )
    crushCandidateVerticalSurfaces(surfaces)
