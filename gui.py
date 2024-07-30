"""
Routines for experimenting with the ideal edge code in Regina's GUI.
"""
from timeit import default_timer
from regina import *
from idealedge import decomposeAlong, idealLoops
from idealedge import isAnnulus, isSphere, fillIdealEdge
from loop import IdealLoop
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


def crushAnnuli( surfaces, threshold=30 ):
    """
    Crushes all annuli in the given packet of normal surfaces, and adds a
    Container of the resulting triangulations as a child of the given packet.

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
        tri.setLabel( "Crushed #{}".format(surfNum) )
        thin = surf.isThinEdgeLink()
        if thin[0] is not None:
            # Adorn label with details of this thin edge link.
            adorn = "Thin edge {}".format( thin[0].index() )
            if thin [1] is not None:
                adorn += " and {}".format( thin[1].index() )
            tri.setLabel( tri.adornedLabel(adorn) )
        components = []
        results.insertChildLast(tri)
        idEdgeDetails = idealLoops(surf)
        if idEdgeDetails:
            # There is only one ideal loop, given by a length-1 sequence of
            # ideal edges.
            idEdge = idEdgeDetails[0][0]
        else:
            idEdge = None
        idComp = None
        if tri.isEmpty():
            tri.setLabel( tri.label() + ": Empty" )
        else:
            if tri.isConnected():
                components.append(tri)
                if idEdge is not None:
                    idComp = 0
            else:
                tri.setLabel( tri.label() + ": Disconnected" )
                for compNum, c in enumerate( tri.triangulateComponents() ):
                    comp = PacketOfTriangulation3(c)
                    comp.setLabel( "Component #{}".format(compNum) )
                    components.append(comp)
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
                comp.setLabel( comp.label() + ": INVALID" )

                # Fill in invalid boundary.
                filled = PacketOfTriangulation3(comp)
                invIdEdge = fillIdealEdge(filled)
                filled.setLabel( comp.adornedLabel(
                    "Closed, ideal edge {}".format(
                        invIdEdge.index() ) ) )
                comp.insertChildLast(filled)

                # Just in case, let's see if we can simplify and identify the
                # manifold given by drilling out the ideal edge.
                drilled = PacketOfTriangulation3(filled)
                filled.insertChildLast(drilled)
                ide = drilled.edge( invIdEdge.index() )
                drilled.setLabel( comp.adornedLabel(
                    "Closed, pinched edge {}".format( ide.index() ) ) )
                drilled.pinchEdge(ide)
                drilled.intelligentSimplify()
                drilled.intelligentSimplify()
                if ( ( drilled.knowsSolidTorus() or
                    drilled.size() < threshold ) and
                    drilled.isSolidTorus() ):
                    name = "Ideal solid torus"
                else:
                    # Try to combinatorially recognise after truncating the
                    # ideal vertex.
                    trunc = PacketOfTriangulation3(drilled)
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
                    trunc.setLabel( drilled.adornedLabel(
                        "Truncated" ) + ": {}".format(name) )
                drilled.setLabel(
                        drilled.label() + ": {}".format(name) )

                # Decompose the filled manifold into prime pieces (unless it
                # has too many tetrahedra).
                print( "        Attempted prime decomposition: {}.".format(
                    recogniseSummands( filled, threshold ) ) )
            else:
                # If this component contains the ideal edge, then attempt to
                # simplify (and possibly identify) the drilled manifold.
                if compNum == idComp:
                    drilled = PacketOfTriangulation3(comp)
                    comp.insertChildLast(drilled)
                    ide = drilled.tetrahedron( idEdge[0] ).edge( idEdge[1] )

                    # Need to label *before* drilling.
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
                        trunc.setLabel( drilled.adornedLabel(
                            "Truncated" ) + ": {}".format(name) )
                    drilled.setLabel(
                            drilled.label() + ": {}".format(name) )

                # Decompose this component into prime pieces (unless this
                # component has too many tetrahedra).
                print( "        Attempted prime decomposition: {}.".format(
                    recogniseSummands( comp, threshold ) ) )

    # All done!
    print()
    print( "Time: {:.6f}. All done!".format(
        default_timer() - start ) )
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
