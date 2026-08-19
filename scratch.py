"""
Scratch work for ideal edges.
"""
from sys import argv
from timeit import default_timer
from regina import *
from idealedge import decomposeAlong, newIdealLoopEmbs, fillIdealEdges
from loop import IdealLoop, BoundsDisc
from pinch import drillMeridian
from wedge import wedgeCycles
from construct.sfs import orientableSFS
from aux.tetrenum import tetRenumbering
from aux.quad import tetHasQuads
from aux.surface import isSphere, isAnnulus
#TODO Replace wedgeCycles() with nonSurvivingTriangularOrbitCounts()


#TODO meridian() is never used anywhere. Just delete it?
def meridian( tri, edgeIndex ):
    """
    Drills out an edge loop e (corresponding to the given triangulation and
    edge index), and returns the resulting TriangulationWithBoundaryLoops,
    which will have a single BoundaryLoop corresponding to the meridian of the
    drilled loop.

    Pre-condition:
    --> The edge given by tri.edge(edgeIndex) must lie entirely in the
        interior of tri, and the two endpoints of this edge must be
        identified.
    """
    return drillMeridian( IdealLoop( [ tri.edge(edgeIndex) ] ) )


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


#TODO Experiment with crushing Mobius bands as well.
def crushAnnuli( surfaces, threshold=30 ):
    """
    Crushes all annuli in the given list of normal surfaces.

    If the given surfaces are contained in a PacketOfNormalSurface, then this
    routine adds a Container of the crushed triangulations as a child of the
    given surfaces packet. Otherwise, this routine simply prints details of
    the crushed triangulations.

    This routine attempts to identify the topology of the manifold that
    results from crushing. The main strategy is to simplify and attempt
    combinatorial recognition. Additionally, whenever this routine encounters
    a component whose number of tetrahedra is strictly less than the
    threshold (default 30), it will also use more computationally intensive
    recognition algorithms involving normal surfaces.

    Pre-condition:
    --> For each annulus in the given list of surfaces, every boundary
        component incident to this annulus must be a two-triangle torus.
    """
    start = default_timer()
    usingPackets = isinstance( surfaces, PacketOfNormalSurfaces )
    if usingPackets:
        results = Container( "Crush annuli" )
        surfaces.insertChildLast(results)
    annulusCount = 0
    for surfNum, surf in enumerate(surfaces):
        if not isAnnulus(surf):
            continue
        thin = surf.isThinEdgeLink()
        if thin[0] is not None:
            # Don't bother with thin edge links.
            continue
        annulusCount += 1
        print()
        print( "Time: {:.6f}. Crush #{}.".format(
            default_timer() - start, surfNum) )

#        # Is the current annulus a thin edge link?
#        thin = surf.isThinEdgeLink()
#        if thin[0] is None:
#            thinAdorn = ""

        # Crush, and find the ideal edge amongst the components of the
        # resulting triangulation.
        crushedName = "Crushed #{}".format(surfNum)
        for _, twist in wedgeCycles(surf):
            if twist != 0:
                crushedName += " (Lost (3,{}))".format(twist)
        #NOTE Crushing preserves orientation.
        tri = PacketOfTriangulation3( surf.crush() )
        try:
            filledH1 = filledHomology(surf)
        except ValueError:
            filledH1 = "unknown"
        print( "--> Dehn-filled homology: {} -> {}".format(
            surf.triangulation().homology(), filledH1 ) )
        if usingPackets:
            tri.setLabel(crushedName)
            results.insertChildLast(tri)
        else:
            print(crushedName)
#        thin = surf.isThinEdgeLink()
#        if thin[0] is not None:
#            # Adorn label with details of this thin edge link.
#            adorn = "Thin edge {}".format( thin[0].index() )
#            if thin [1] is not None:
#                adorn += " and {}".format( thin[1].index() )
#            if usingPackets:
#                tri.setLabel( tri.adornedLabel(adorn) )
#            else:
#                # Or just print if we're not using packets.
#                print(adorn)
        components = []
        idEdgeDetailsInOldTri = newIdealLoopEmbs(surf)
        if idEdgeDetailsInOldTri:
            # There is only one ideal loop, given by a length-1 sequence of
            # ideal edges.
            oldEmb = idEdgeDetailsInOldTri[0][0]

            # Translate oldEmb into an edge embedding in the crushed
            # triangulation tri.
            doomed = [ tet for tet in surf.triangulation().tetrahedra()
                      if tetHasQuads( surf, tet.index() ) ]
            tetIndicesAfterCrush = tetRenumbering(doomed)
            crushedTet = tri.tetrahedron(
                    tetIndicesAfterCrush[ oldEmb.tetrahedron().index() ] )
            idEdgeEmb = EdgeEmbedding3( crushedTet, oldEmb.vertices() )
        else:
            idEdgeEmb = None

        # Split tri into components, and find the ideal edge amongst the newly
        # split components.
        idComp = None
        if tri.isEmpty():
            if usingPackets:
                tri.setLabel( tri.label() + ": Empty" )
            else:
                print("Empty triangulation")
        else:
            if tri.isConnected():
                components.append(tri)
                if idEdgeEmb is not None:
                    idComp = 0
            else:
                if usingPackets:
                    tri.setLabel( tri.label() + ": Disconnected" )
                else:
                    print("Disconnected triangulation")
                #NOTE Splitting into components naturally preserves orientation.
                for compNum, c in enumerate( tri.triangulateComponents() ):
                    comp = PacketOfTriangulation3(c)
                    components.append(comp)
                    if usingPackets:
                        comp.setLabel( "Component #{}".format(compNum) )
                        tri.insertChildLast(comp)

                # Find the component containing the ideal edge, and adjust
                # the ideal tetrahedron index.
                if idEdgeEmb is not None:
                    idComp = idEdgeEmb.tetrahedron().component().index()
                    idTeti = 0
                    for tet in tri.tetrahedra():
                        if tet.component().index() == idComp:
                            if tet == idEdgeEmb.tetrahedron():
                                idEdgeInfo = ( idTeti, idEdgeEmb.vertices() )
                                break
                            else:
                                idTeti += 1

        # Go through the components and try to identify their topology.
        for compNum, comp in enumerate(components):
            print( "    Time: {:.6f}. Component #{}.".format(
                default_timer() - start, compNum ) )
            if not comp.isValid():
                if usingPackets:
                    comp.setLabel( comp.label() + ": INVALID" )
                else:
                    print( "        INVALID" )

                # Fill in invalid boundary.
                #NOTE The mess below is to work around layerOn() failing with
                #   PacketOfTriangulation3.
                filled = Triangulation3(comp)
                #NOTE fillIdealEdges() preserves orientation.
                endpoints = set()
                for v in filled.vertices():
                    if not v.isValid():
                        endpoints.add( v.index() )
                invIdEdgeIndex = fillIdealEdges(
                        filled, endpoints )[0].index()
                filled = PacketOfTriangulation3(filled)
                print( "--> filled component homology: {}".format(
                    filled.homology() ) )
                invIdEdge = filled.edge(invIdEdgeIndex)
                if usingPackets:
                    filled.setLabel( comp.adornedLabel(
                        "Closed, ideal edge {}".format(
                            invIdEdge.index() ) ) )
                    comp.insertChildLast(filled)
                else:
                    print( "        Closed, ideal edge {}".format(
                        invIdEdge.index() ) )

                # Have we isolated a single exceptional fibre?
                invIdLoop = IdealLoop( [invIdEdge] )
                try:
                    # The meridian of the ideal loop is a candidate for a
                    # regular fibre.
                    #NOTE Drilling preserves orientation.
                    triWithMeridian = drillMeridian(invIdLoop)
                except BoundsDisc:
                    # The meridian bounds a disc "on the outside", so the
                    # filled triangulation must have been S2 x S1. In
                    # particular, the meridian cannot be a regular fibre.
                    if usingPackets:
                        filled.setLabel(
                                filled.label() + ": {}".format(
                                    "S2 x S1, meridian is not a fibre" ) )
                    else:
                        print( "        S2 x S1, meridian is not a fibre" )
                else:
                    # Successfully drilled.
                    #NOTE Simplification preserves orientation.
                    triWithMeridian.minimiseBoundary()
                    triWithMeridian.simplify()
                    triWithMeridian.simplify()
                    drilled = PacketOfTriangulation3(
                            triWithMeridian.triangulation() )
                    print( "Oriented? {}".format( drilled.isOriented() ) )
                    if usingPackets:
                        filled.insertChildLast(drilled)

                    # There is only one BoundaryLoop, corresponding to the
                    # meridian. Also, because we minimised the boundary, the
                    # meridian is guaranteed to be given by a single edge.
                    merEdgeIndex = triWithMeridian[0][0]
                    if usingPackets:
                        drilled.setLabel( comp.adornedLabel(
                            "Drilled, meridian edge {}".format(merEdgeIndex) ) )
                    else:
                        print( "        Drilled, meridian edge {} (Time: {:.6f})".format(
                            merEdgeIndex, default_timer() - start ) )

                    # If the drilled triangulation is a solid torus, then
                    # finding the compression disc D will tell us the
                    # parameters of the exceptional fibre.
                    #
                    # In detail, let M denote the weight of D on the
                    # meridian and let E denote the weight of D on one of the
                    # other boundary edges (labelled e in the diagram below).
                    # Orient the meridian edge (upwards in the diagram below)
                    # and number the intersection points in order from 0 to
                    # M-1. An arc of the boundary of D leaving point p along
                    # the meridian will return to the meridian at:
                    #       (p plus/minus E) mod M
                    # The choice between p+E or p-E depends on the direction
                    # of the arc, as well as on whether E > M or M > E.
                    #
                    #           e
                    #       +-------+
                    #       |       |
                    #   mer ^       ^
                    #       |       |
                    #       +-------+
                    #
                    # Thus, ignoring orientation, we can determine the
                    # parameters of the exceptional fibre by computing the
                    # multiplicative inverse of E mod M (which exists because
                    # gcd(E,M) = 1).
                    surf = drilled.nonTrivialSphereOrDisc()
                    if surf is None:
                        # No compression disc means we have not yet cut out a
                        # single fibre.
                        # Continue by decomposing along spheres.
                        name = "Decomposed into fibres: "
                        decomposedLoops = decomposeAlongSpheres(invIdLoop)
                        for newLoop in decomposedLoops:
                            # We should be able to drill and get Seifert
                            # fibre parameters.
                            try:
                                newTriWithMeridian = drillMeridian(newLoop)
                            except BoundsDisc:
                                name += "N/A (S2 x S1); "
                            else:
                                newTriWithMeridian.minimiseBoundary()
                                newTriWithMeridian.simplify()
                                newTriWithMeridian.simplify()
                                newDrilled = PacketOfTriangulation3(
                                        newTriWithMeridian.triangulation() )

                                # There is only one BoundaryLoop, corresponding to the
                                # meridian. Also, because we minimised the boundary, the
                                # meridian is guaranteed to be given by a single edge.
                                newMerEdgeIndex = newTriWithMeridian[0][0]
                                newSurf = newDrilled.nonTrivialSphereOrDisc()
                                if newSurf is None:
                                    name += "N/A (no disc); "
                                elif newSurf.eulerChar() == 2:
                                    name += "unknown (found sphere); "
                                else:
                                    name += "({},{}); ".format(
                                            *fibreParams( newSurf, newMerEdgeIndex ) )

                        # Format name correctly.
                        name = name[:-2]
                    elif surf.eulerChar() == 2:
                        #TODO Sphere. Probably want to crush.
                        name = "Contains nontrivial sphere"
                    else:
                        # Use boundary edge weights of the disc to calculate
                        # Seifert parameters.
                        name = "Seifert fibre (p,q)=({},{})".format(
                                *fibreParams( surf, merEdgeIndex ) )
                    if usingPackets:
                        drilled.setLabel(
                                drilled.label() + ": {}".format(name) )
                    else:
                        print( "        " + name )

                #TODO If we didn't get down to a single fibre, then we should
                #   continue by decomposing along spheres that intersect the
                #   ideal loop twice.

#                # Just in case, let's see if we can simplify and identify the
#                # manifold given by drilling out the ideal edge.
#                drilled = PacketOfTriangulation3(filled)
#                filled.insertChildLast(drilled)
#                ide = drilled.edge( invIdEdge.index() )
#                drilled.setLabel( comp.adornedLabel(
#                    "Closed, pinched edge {}".format( ide.index() ) ) )
#                drilled.pinchEdge(ide)
#                drilled.simplify()
#                drilled.simplify()
#                if ( ( drilled.knowsSolidTorus() or
#                    drilled.size() < threshold ) and
#                    drilled.isSolidTorus() ):
#                    name = "Ideal solid torus"
#                else:
#                    # Try to combinatorially recognise after truncating the
#                    # ideal vertex.
#                    trunc = PacketOfTriangulation3(drilled)
#                    drilled.insertChildLast(trunc)
#                    trunc.idealToFinite()
#                    trunc.simplify()
#                    trunc.simplify()
#                    std = StandardTriangulation.recognise(trunc)
#                    if std is None:
#                        name = "Not recognised"
#                        if drilled.knowsSolidTorus():
#                            name += ", not solid torus"
#                    else:
#                        name = std.manifold().name()
#                    trunc.setLabel( drilled.adornedLabel(
#                        "Truncated" ) + ": {}".format(name) )
#                drilled.setLabel(
#                        drilled.label() + ": {}".format(name) )
#
#                # Decompose the filled manifold into prime pieces (unless it
#                # has too many tetrahedra).
#                print( "        Attempted prime decomposition: {}.".format(
#                    recogniseSummands( filled, threshold ) ) )
            else:
                #TODO Experiment with drillMeridian() instead of pinchEdge().
                # If this component contains the ideal edge, then attempt to
                # simplify (and possibly identify) the drilled manifold.
                print( "--> component homology: {}".format(
                    comp.homology() ) )
                if compNum == idComp:
                    idTeti, idVer = idEdgeInfo
                    ide = comp.tetrahedron(idTeti).edge(
                            idVer[0], idVer[1] )
                    idLoop = IdealLoop( [ide] )
                    if usingPackets:
                        comp.setLabel( comp.adornedLabel(
                            "Ideal edge {}".format( ide.index() ) ) )
                    else:
                        print( "        Ideal edge {}".format(
                            ide.index() ) )
                    try:
                        # The meridian of the ideal loop is a candidate for an
                        # exceptional fibre.
                        #NOTE Drilling preserves orientation.
                        triWithMeridian = drillMeridian(idLoop)
                    except BoundsDisc:
                        # The meridian bounds a disc "on the outside", so the
                        # filled triangulation must have been S2 x S1. In
                        # particular, the meridian cannot be an exceptional
                        # fibre.
                        if usingPackets:
                            comp.setLabel(
                                    comp.label() + ": {}".format(
                                        "S2 x S1, meridian is not a fibre" ) )
                        else:
                            print( "        S2 x S1, meridian is not a fibre" )
                    else:
                        # Successfully drilled.
                        #NOTE Simplification preserves orientation.
                        triWithMeridian.minimiseBoundary()
                        triWithMeridian.simplify()
                        triWithMeridian.simplify()
                        drilled = PacketOfTriangulation3(
                                triWithMeridian.triangulation() )
                        print( "Oriented? {}".format( drilled.isOriented() ) )
                        if usingPackets:
                            comp.insertChildLast(drilled)

                        # There is only one BoundaryLoop, corresponding to the
                        # meridian. Also, because we minimised the boundary, the
                        # meridian is guaranteed to be given by a single edge.
                        merEdgeIndex = triWithMeridian[0][0]
                        if usingPackets:
                            drilled.setLabel( comp.adornedLabel(
                                "Drilled, meridian edge {}".format(merEdgeIndex) ) )
                        else:
                            print( "        Drilled, meridian edge {} (Time: {:.6f})".format(
                                merEdgeIndex, default_timer() - start ) )

                        # If the drilled triangulation is a solid torus, then
                        # finding the compression disc D will tell us the
                        # parameters of the exceptional fibre.
                        #
                        # In detail, let M denote the weight of D on the
                        # meridian and let E denote the weight of D on one of the
                        # other boundary edges (labelled e in the diagram below).
                        # Orient the meridian edge (upwards in the diagram below)
                        # and number the intersection points in order from 0 to
                        # M-1. An arc of the boundary of D leaving point p along
                        # the meridian will return to the meridian at:
                        #       (p plus/minus E) mod M
                        # The choice between p+E or p-E depends on the direction
                        # of the arc, as well as on whether E > M or M > E.
                        #
                        #           e
                        #       +-------+
                        #       |       |
                        #   mer ^       ^
                        #       |       |
                        #       +-------+
                        #
                        # Thus, ignoring orientation, we can determine the
                        # parameters of the exceptional fibre by computing the
                        # multiplicative inverse of E mod M (which exists because
                        # gcd(E,M) = 1).
                        surf = drilled.nonTrivialSphereOrDisc()
                        if surf is None:
                            # No compression disc means we have not yet cut out a
                            # single fibre.
                            # Continue by decomposing along spheres.
                            name = "Decomposed into fibres: "
                            decomposedLoops = decomposeAlongSpheres(idLoop)
                            for newLoop in decomposedLoops:
                                # We should be able to drill and get Seifert
                                # fibre parameters.
                                try:
                                    newTriWithMeridian = drillMeridian(newLoop)
                                except BoundsDisc:
                                    name += "N/A (S2 x S1); "
                                else:
                                    newTriWithMeridian.minimiseBoundary()
                                    newTriWithMeridian.simplify()
                                    newTriWithMeridian.simplify()
                                    newDrilled = PacketOfTriangulation3(
                                            newTriWithMeridian.triangulation() )

                                    # There is only one BoundaryLoop, corresponding to the
                                    # meridian. Also, because we minimised the boundary, the
                                    # meridian is guaranteed to be given by a single edge.
                                    newMerEdgeIndex = newTriWithMeridian[0][0]
                                    newSurf = newDrilled.nonTrivialSphereOrDisc()
                                    if newSurf is None:
                                        name += "N/A (no disc); "
                                    elif newSurf.eulerChar() == 2:
                                        name += "unknown (found sphere); "
                                    else:
                                        name += "({},{}); ".format(
                                                *fibreParams( newSurf, newMerEdgeIndex ) )

                            # Format name correctly.
                            name = name[:-2]
                        elif surf.eulerChar() == 2:
                            #TODO Sphere. Probably want to crush.
                            name = "Contains nontrivial sphere"
                        else:
                            # Use boundary edge weights of the disc to calculate
                            # Seifert parameters.
                            name = "Seifert fibre (p,q)=({},{})".format(
                                    *fibreParams( surf, merEdgeIndex ) )
                        if usingPackets:
                            drilled.setLabel(
                                    drilled.label() + ": {}".format(name) )
                        else:
                            print( "        " + name )

                    #TODO If we didn't get down to a single fibre, then we should
                    #   continue by decomposing along spheres that intersect the
                    #   ideal loop twice.

#                    drilled = PacketOfTriangulation3(comp)
#                    if usingPackets:
#                        comp.insertChildLast(drilled)
#                    ide = drilled.tetrahedron( idEdge[0] ).edge( idEdge[1] )
#
#                    # Need to label *before* drilling.
#                    if usingPackets:
#                        drilled.setLabel( comp.adornedLabel(
#                            "Pinched edge {}".format( ide.index() ) ) )
#                        comp.setLabel( comp.adornedLabel(
#                            "Ideal edge {}".format( ide.index() ) ) )
#                    drilled.pinchEdge(ide)
#                    drilled.simplify()
#                    drilled.simplify()
#
#                    # Try to recognise the drilled manifold.
#                    if ( ( drilled.knowsSolidTorus() or
#                        drilled.size() < threshold ) and
#                        drilled.isSolidTorus() ):
#                        name = "Ideal solid torus"
#                    else:
#                        # Try to combinatorially recognise after truncating
#                        # the ideal vertex.
#                        trunc = PacketOfTriangulation3(drilled)
#                        if usingPackets:
#                            drilled.insertChildLast(trunc)
#                        trunc.idealToFinite()
#                        trunc.simplify()
#                        trunc.simplify()
#                        std = StandardTriangulation.recognise(trunc)
#                        if std is None:
#                            name = "Not recognised"
#                            if drilled.knowsSolidTorus():
#                                name += ", not solid torus"
#                        else:
#                            name = std.manifold().name()
#                        if usingPackets:
#                            trunc.setLabel( drilled.adornedLabel(
#                                "Truncated" ) + ": {}".format(name) )
#                    if usingPackets:
#                        drilled.setLabel(
#                                drilled.label() + ": {}".format(name) )
#                    else:
#                        print( "        " + name)
#
#                # Decompose this component into prime pieces (unless this
#                # component has too many tetrahedra).
#                print( "        Attempted prime decomposition: {}.".format(
#                    recogniseSummands( comp, threshold ) ) )

    # All done!
    print()
    print( "Time: {:.6f}. All done!".format(
        default_timer() - start ) )
    if usingPackets:
        results.setLabel( results.adornedLabel(
            "Total {}".format(annulusCount) ) )


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


def crushSpheres( surfaces, idealEdgeIndex, threshold=30 ):
    """
    """
    results = Container( "Decompose along 2-spheres" )
    surfaces.insertChildLast(results)
    for surfNum, surf in enumerate(surfaces):
        if not isSphere(surf):
            continue
        try:
            #TODO This needs to be updated.
            #TODO Update usage to account for extra book-keeping.
            pieces = decomposeAlong( surf, {idealEdgeIndex} )
        except ValueError:
            continue
        container = Container( "Decompose along #{}".format(surfNum) )
        results.insertChildLast(container)
        for i, piece in enumerate(pieces):
            tri = PacketOfTriangulation3( piece[0] )
            loops = piece[1]
            tri.setLabel( "Component #{}: {}".format(
                i, loops ) )
            container.insertChildLast(tri)

            # Is tri a 3-sphere?
            if ( tri.knowsSphere() or tri.size() < threshold ):
                if tri.isSphere():
                    name = "S3"
                else:
                    name = "Not S3"
            else:
                name = "Not recognised"
            tri.setLabel( tri.label() + ": {}".format(name) )

            # Build drilled 3-manifold.
            drilled = PacketOfTriangulation3(tri)
            drilled.setLabel( tri.adornedLabel(
                "Pinched ideal edges" ) )
            tri.insertChildLast(drilled)
            for t, e in loops:
                drilled.pinchEdge( drilled.tetrahedron(t).edge(e) )
                drilled.simplify()
                drilled.simplify()

                # Is drilled a solid torus?
                if ( drilled.knowsSolidTorus() or
                        drilled.size() < threshold ):
                    if drilled.isSolidTorus():
                        name = "Ideal solid torus"
                    else:
                        name = "Ideal, not solid torus"
                else:
                    name = "Ideal, not recognised"
                drilled.setLabel(
                        drilled.label() + ": {}".format(name) )
            #TODO
            pass
        #TODO
        pass
    #TODO
    return


def decomposeAlongSpheres(idealLoop):
    """
    Returns a list of ideal loops given by repeatedly decomposing the given
    ideal loop along spheres that intersect the ideal loop twice.
    """
    # Repeatedly crush along spheres that intersect the ideal loop at most
    # twice.
    toProcess = [idealLoop]
    decomposedLoops = []
    while toProcess:
        oldLoop = toProcess.pop()
        tri = oldLoop.triangulation()

        # Search for a suitable sphere to crush.
        enumeration = TreeEnumeration( tri, NS_QUAD )
        while True:
            # Get the next 2-sphere.
            if enumeration.next():
                sphere = enumeration.buildSurface()
                if not isSphere(sphere):
                    continue
            else:
                # No suitable 2-sphere means that we're done with the current
                # oldLoop.
                decomposedLoops.append(oldLoop)
                break

            # Does the sphere intersect the ideal loop at most twice?
            wt = oldLoop.weight(sphere)
            if wt != 2:
                #TODO Actually do something with the following cases.
                if wt == 0:
                    print( "Found sphere disjoint from ideal loop!" )
                if wt == 1:
                    print( "Found sphere intersecting ideal loop once!" )

                # Continue searching for suitable spheres.
                continue

            # See what happens if we crush.
            lostFibres = ""
            for _, twist in wedgeCycles(sphere):
                if twist != 0:
                    lostFibres += " Lost (3,{}).".format(twist)
            if lostFibres:
                print( lostFibres[1:] )
            #TODO Update usage to account for extra book-keeping.
            decomposed = decomposeAlong( sphere, [oldLoop] )
            for newLoops in decomposed:
                if newLoops:
                    # We are guaranteed to have len(newLoops) == 1.
                    toProcess.append( newLoops[0] )

            # Found and crushed a suitable sphere, so stop enumerating.
            break

    # If we reach this point, then we have decomposed as far as possible, and
    # everything remaining has no suitable spheres.
    return decomposedLoops


def recogniseSummands( tri, threshold=40 ):
    """
    Attempts to recognise the prime summands of the given triangulation.

    This routine only proceeds with performing the prime decomposition if the
    number of tetrahedra in tri is strictly less than the threshold (default
    40), and returns True if and only if this is the case.
    """
    if tri.size() >= threshold:
        return False
    summands = tri.summands()
    if len(summands) == 0:
        tri.setLabel( tri.label() + ": S3" )
    elif len(summands) == 1:
        # Try combinatorial recognition.
        std = StandardTriangulation.recognise( summands[0] )
        if std is None:
            name = "Prime, not recognised"
        else:
            name = std.manifold().name()
        tri.setLabel( tri.label() + ": {}".format(name) )
    else:
        tri.setLabel( tri.label() + ": Non-prime" )

        # Find *all* quad vertex normal 2-spheres.
        surfs = NormalSurfaces( tri, NS_QUAD, NS_VERTEX )
        sphereFilter = SurfaceFilterProperties()
        sphereFilter.setEulerChars( [2] )
        sphereFilter.setCompactness( BoolSet(True) )
        sphereFilter.setOrientability( BoolSet(True) )
        sphereFilter.setRealBoundary( BoolSet(False) )
        spheres = PacketOfNormalSurfaces( surfs, sphereFilter )
        spheres.setLabel( "Quad vertex 2-spheres (Total: {})".format(
            spheres.size() ) )
        tri.insertChildLast(spheres)

        # Classify the summands.
        sumContainer = Container( "Summands (Total: {})".format(
            len(summands) ) )
        tri.insertChildLast(sumContainer)
        for sumNum, s in enumerate(summands):
            summand = PacketOfTriangulation3(s)
            sumContainer.insertChildLast(summand)

            # Try to combinatorially recognise this summand.
            std = StandardTriangulation.recognise(summand)
            if std is None:
                name = "Prime, not recognised"
            else:
                name = std.manifold().name()
            summand.setLabel( "Summand #{}: {}".format(
                sumNum, name ) )
    return True


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
#    params = [ int(n) for n in argv[1:] ]
#    manifold = SFSpace()
#    manifold.insertFibre(3,1)
#    while params:
#        q = params.pop()
#        p = params.pop()
#        manifold.insertFibre(p,q)
#    tri = manifold.construct()
#    tri.removeTetrahedronAt(3)
#    tri.orient()
#    tri.simplify()
#    tri.simplify()
##    p = int( argv[1] )
##    q = int( argv[2] )
##    knot = ExampleLink.torus(p,q)
##    ext = knot.complement()
##    ext.idealToFinite()
##    ext.simplify()
##    ext.simplify()
##    surfaces = NormalSurfaces( ext, NS_QUAD, NS_VERTEX )
    #NOTE As of Regina 7.4, NS_QUAD and NS_VERTEX have been deprecated and
    #       replaced with NormalCoords.Quad and NormalList.Vertex,
    #       respectively.
    surfaces = NormalSurfaces( tri, NormalCoords.Quad, NormalList.Vertex )
    crushAnnuli(surfaces)
