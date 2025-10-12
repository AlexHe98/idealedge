"""
Routines for experimenting with the ideal edge code in Regina's GUI.
"""
from sys import argv
from timeit import default_timer
from regina import *
from idealedge import decomposeAlong, idealLoops
from idealedge import isAnnulus, isSphere, fillIdealEdge
from loop import IdealLoop, BoundsDisc
from subdivide import drillMeridian


def meridian( tri, edgeIndex ):
    """
    Drills out an edge loop e (corresponding to the given triangulation and
    edge index), and returns the resulting meridian curve.

    Pre-condition:
    --> The edge given by tri.edge(edgeIndex) must lie entirely in the
        interior of tri, and the two endpoints of this edge must be
        identified.
    """
    return drillMeridian( IdealLoop( [ tri.edge(edgeIndex) ] ) )


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
        annulusCount += 1
        print()
        print( "Time: {:.6f}. Crush #{}.".format(
            default_timer() - start, surfNum) )

        # Is the current annulus a thin edge link?
        thin = surf.isThinEdgeLink()
        if thin[0] is None:
            thinAdorn = ""

        # Crush, and find the ideal edge amongst the components of the
        # resulting triangulation.
        tri = PacketOfTriangulation3( surf.crush() )
        if usingPackets:
            tri.setLabel( "Crushed #{}".format(surfNum) )
            results.insertChildLast(tri)
        thin = surf.isThinEdgeLink()
        if thin[0] is not None:
            # Adorn label with details of this thin edge link.
            adorn = "Thin edge {}".format( thin[0].index() )
            if thin [1] is not None:
                adorn += " and {}".format( thin[1].index() )
            if usingPackets:
                tri.setLabel( tri.adornedLabel(adorn) )
            else:
                # Or just print if we're not using packets.
                print(adorn)
        components = []
        idEdgeDetails = idealLoops(surf)
        if idEdgeDetails:
            # There is only one ideal loop, given by a length-1 sequence of
            # ideal edges.
            idEdge = idEdgeDetails[0][0]
        else:
            idEdge = None
        idComp = None
        if tri.isEmpty():
            if usingPackets:
                tri.setLabel( tri.label() + ": Empty" )
            else:
                print("Empty triangulation")
        else:
            if tri.isConnected():
                components.append(tri)
                if idEdge is not None:
                    idComp = 0
            else:
                if usingPackets:
                    tri.setLabel( tri.label() + ": Disconnected" )
                else:
                    print("Disconnected triangulation")
                for compNum, c in enumerate( tri.triangulateComponents() ):
                    comp = PacketOfTriangulation3(c)
                    components.append(comp)
                    if usingPackets:
                        comp.setLabel( "Component #{}".format(compNum) )
                        tri.insertChildLast(comp)

                # Find the component containing the ideal edge, and adjust
                # the ideal tetrahedron index.
                if idEdge is not None:
                    idComp = tri.tetrahedron( idEdge[0] ).component().index()
                    idTeti = 0
                    for tet in tri.tetrahedra():
                        if tet.component().index() == idComp:
                            if tet.index() == idEdge[0]:
                                idEdge = ( idTeti, idEdge[1] )
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
                filled = PacketOfTriangulation3(comp)
                invIdEdge = fillIdealEdge(filled)
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
                    # The meridian of the ideal loop is a candidate for an
                    # exceptional fibre.
                    mer = drillMeridian(invIdLoop)
                except BoundsDisc:
                    # The meridian bounds a disc "on the outside", so the
                    # filled triangulation must have been S2 x S1. In
                    # particular, the meridian cannot be an exceptional
                    # fibre.
                    if usingPackets:
                        filled.setLabel(
                                filled.label() + ": {}".format(
                                    "S2 x S1, meridian is not a fibre" ) )
                    else:
                        print( "        S2 x S1, meridian is not a fibre" )
                else:
                    # Successfully drilled.
                    mer.minimiseBoundary()
                    mer.simplify()
                    mer.simplify()
                    drilled = PacketOfTriangulation3( mer.triangulation() )
                    if usingPackets:
                        filled.insertChildLast(drilled)

                    # Because we minimised the boundary, the meridian is
                    # guaranteed to be given by a single edge.
                    merEdgeIndex = mer[0]
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
                        name = "Not a fibred solid torus"
                    elif surf.eulerChar() == 2:
                        #TODO Sphere. Probably want to crush.
                        name = "Contains nontrivial sphere"
                    else:
                        # Use boundary edge weights of the disc to calculate
                        # Seifert parameters (as outlined above).
                        merWt = surf.edgeWeight(merEdgeIndex).safeLongValue()
                        for e in drilled.edges():
                            if e.index() == merEdgeIndex or not e.isBoundary():
                                continue

                            # Found another boundary edge.
                            bdyWt = surf.edgeWeight( e.index() ).safeLongValue()
                            break
                        name = "Seifert fibre (p,q)=({},{})".format(
                                merWt, pow( bdyWt, -1, merWt ) )
                    if usingPackets:
                        drilled.setLabel(
                                drilled.label() + ": {}".format(name) )
                    else:
                        print( "        " + name )
                #TODO

#                # Just in case, let's see if we can simplify and identify the
#                # manifold given by drilling out the ideal edge.
#                drilled = PacketOfTriangulation3(filled)
#                filled.insertChildLast(drilled)
#                ide = drilled.edge( invIdEdge.index() )
#                drilled.setLabel( comp.adornedLabel(
#                    "Closed, pinched edge {}".format( ide.index() ) ) )
#                drilled.pinchEdge(ide)
#                drilled.intelligentSimplify()
#                drilled.intelligentSimplify()
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
#                    trunc.intelligentSimplify()
#                    trunc.intelligentSimplify()
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
                if compNum == idComp:
                    drilled = PacketOfTriangulation3(comp)
                    if usingPackets:
                        comp.insertChildLast(drilled)
                    ide = drilled.tetrahedron( idEdge[0] ).edge( idEdge[1] )

                    # Need to label *before* drilling.
                    if usingPackets:
                        drilled.setLabel( comp.adornedLabel(
                            "Pinched edge {}".format( ide.index() ) ) )
                        comp.setLabel( comp.adornedLabel(
                            "Ideal edge {}".format( ide.index() ) ) )
                    drilled.pinchEdge(ide)
                    drilled.intelligentSimplify()
                    drilled.intelligentSimplify()

                    # Try to recognise the drilled manifold.
                    if ( ( drilled.knowsSolidTorus() or
                        drilled.size() < threshold ) and
                        drilled.isSolidTorus() ):
                        name = "Ideal solid torus"
                    else:
                        # Try to combinatorially recognise after truncating
                        # the ideal vertex.
                        trunc = PacketOfTriangulation3(drilled)
                        if usingPackets:
                            drilled.insertChildLast(trunc)
                        trunc.idealToFinite()
                        trunc.intelligentSimplify()
                        trunc.intelligentSimplify()
                        std = StandardTriangulation.recognise(trunc)
                        if std is None:
                            name = "Not recognised"
                            if drilled.knowsSolidTorus():
                                name += ", not solid torus"
                        else:
                            name = std.manifold().name()
                        if usingPackets:
                            trunc.setLabel( drilled.adornedLabel(
                                "Truncated" ) + ": {}".format(name) )
                    if usingPackets:
                        drilled.setLabel(
                                drilled.label() + ": {}".format(name) )
                    else:
                        print( "        " + name)

                # Decompose this component into prime pieces (unless this
                # component has too many tetrahedra).
                print( "        Attempted prime decomposition: {}.".format(
                    recogniseSummands( comp, threshold ) ) )

    # All done!
    print()
    print( "Time: {:.6f}. All done!".format(
        default_timer() - start ) )
    if usingPackets:
        results.setLabel( results.adornedLabel(
            "Total {}".format(annulusCount) ) )


def decomposeAlongSpheres( surfaces, idealEdgeIndex, threshold=30 ):
    """
    """
    results = Container( "Decompose along 2-spheres" )
    surfaces.insertChildLast(results)
    for surfNum, surf in enumerate(surfaces):
        if not isSphere(surf):
            continue
        try:
            #TODO This needs to be updated.
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
                drilled.intelligentSimplify()
                drilled.intelligentSimplify()

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
    p = int( argv[1] )
    q = int( argv[2] )
    knot = ExampleLink.torus(p,q)
    ext = knot.complement()
    ext.idealToFinite()
    ext.intelligentSimplify()
    ext.intelligentSimplify()
    surfaces = NormalSurfaces( ext, NS_QUAD, NS_VERTEX )
    crushAnnuli(surfaces)
