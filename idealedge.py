"""
Find the ideal edge after crushing a separating annulus in a one-vertex
triangulation of a 3-manifold with torus boundary.
"""
from regina import *


def isAnnulus(s):
    """
    Is the given normal surface s an annulus?
    """
    return ( s.isCompact() and s.isOrientable() and
            s.hasRealBoundary() and s.eulerChar() == 0 )


def idealEdge(annulus):
    """
    Returns details of the ideal edge after crushing the given annulus.

    Specifically, if the ideal edge belongs to a component that gets
    destroyed after crushing, then this routine returns None. Otherwise, this
    routine returns a pair (t,e), where t is the index (after crushing) of a
    tetrahedron that will meet the ideal edge after crushing, and e is an
    edge number of this tetrahedron that corresponds to the ideal edge.

    Warning:
    --> This routine assumes that the annulus is a separating surface in a
        one-vertex triangulation of a 3-manifold with torus boundary. This
        routine does not check that these conditions are satisfied.
    """
    tri = annulus.triangulation()

    # Use annulus to divide edges of tri into segments.
    segments = []
    for e in tri.edges():
        ei = e.index()
        wt = annulus.edgeWeight(ei).safeLongValue()
        for s in range( wt + 1 ):
            segments.append( (ei,s) )

    # Find "target" segments.
    targets = dict()
    shift = 0
    for tet in tri.tetrahedra():
        teti = tet.index()
        hasQuads = False
        for q in range(3):
            if annulus.quads( teti, q ).safeLongValue() > 0:
                hasQuads = True
                break
        if hasQuads:
            shift += 1
            continue

        # No quads in tet, so there is a cell in the centre that survives
        # crushing. Find the edges of this cell that survive.
        for en in range(6):
            v = tet.edgeMapping(en)[0]
            ei = tet.edge(en).index()
            s = annulus.triangles( teti, v ).safeLongValue()
            targets[ (ei,s) ] = ( teti - shift, en )

    # Starting segment for depth-first search.
    for e in tri.edges():
        ei = e.index()
        if ( e.isBoundary() and
                annulus.edgeWeight(ei).safeLongValue() >= 2 ):
            start = ( ei, 1 )

            # If the starting segment is one of the targets, then we are
            # already done.
            output = targets.get( start, None )
            if output is not None:
                return output

    # Find ideal edge using depth-first search.
    stack = [start]
    visited = set()
    while stack:
        current = stack.pop()
        if current in visited:
            continue

        # We haven't visited the current segment yet, so we need to find all
        # segments that are adjacent to it along parallel cells.
        ei, seg = current
        e = tri.edge(ei)
        wt = annulus.edgeWeight(ei).safeLongValue()
        visited.add(current)    # Record as visited now, so we don't forget.
        for emb in e.embeddings():
            tet = emb.tetrahedron()
            teti = tet.index()
            en = emb.face()
            ver = emb.vertices()

            # To locate the relevant parallel cells in tet, need to get the
            # normal coordinates incident to e.
            f = [ annulus.triangles( teti, ver[i] ).safeLongValue()
                    for i in range(2) ]
            q = 0
            qType = None
            for qt in range(3):
                if qt in { en, 5-en }:
                    continue
                quads = annulus.quads( teti, qt ).safeLongValue()
                if quads > 0:
                    q = quads
                    qType = qt
                    break

            # Does the current segment belong to a parallel cell in this tet?
            # If so, then we need to find all adjacent segments.
            if seg == f[0] or seg == f[0] + q:
                continue
            if seg < f[0]:
                # The parallel cell is a triangular cell at vertex ver[0].
                for enOther in range(6):
                    if enOther == en:
                        continue
                    eiOther = tet.edge(enOther).index()
                    verOther = tet.edgeMapping(enOther)

                    # Is the current segment adjacent to a segment of the
                    # edge numbered enOther?
                    if verOther[0] == ver[0]:
                        adjacent = ( eiOther, seg )
                    elif verOther[1] == ver[0]:
                        wtOther = annulus.edgeWeight(eiOther).safeLongValue()
                        adjacent = ( eiOther, wtOther - seg )
                    else:
                        continue

                    # If the adjacent segment is one of the targets, then we
                    # are done; otherwise, we add it to the stack.
                    output = targets.get( adjacent, None )
                    if output is not None:
                        return output
                    else:
                        stack.append(adjacent)
            elif q > 0 and seg < f[0] + q:
                # The parallel cell is a quadrilateral cell.
                # This quadrilateral cell divides tet into two "sides".
                side0 = { 0, qType + 1 }
                if ver[0] not in side0:
                    side0 = {0,1,2,3} - side0
                qDepth = seg - f[0]     # 1 <= qDepth <= q - 1
                for enOther in range(6):
                    if enOther in { en, qType, 5 - qType }:
                        continue
                    eiOther = tet.edge(enOther).index()
                    verOther = tet.edgeMapping(enOther)

                    # The current segment is adjacent to a segment of the
                    # edge numbered enOther. Find this adjacent segment.
                    fOther = annulus.triangles(
                            teti, verOther[0] ).safeLongValue()
                    if verOther[0] in side0:
                        adjacent = ( eiOther, fOther + qDepth )
                    else:
                        adjacent = ( eiOther, fOther + q - qDepth )

                    # If the adjacent segment is one of the targets, then we
                    # are done; otherwise, we add it to the stack.
                    output = targets.get( adjacent, None )
                    if output is not None:
                        return output
                    else:
                        stack.append(adjacent)
            else:
                # The parallel cell is a triangular cell at vertex ver[1].
                for enOther in range(6):
                    if enOther == en:
                        continue
                    eiOther = tet.edge(enOther).index()
                    verOther = tet.edgeMapping(enOther)

                    # Is the current segment adjacent to a segment of the
                    # edge numbered enOther?
                    if verOther[0] == ver[1]:
                        adjacent = ( eiOther, wt - seg )
                    elif verOther[1] == ver[1]:
                        wtOther = annulus.edgeWeight(eiOther).safeLongValue()
                        adjacent = ( eiOther, wtOther - wt + seg )
                    else:
                        continue

                    # If the adjacent segment is one of the targets, then we
                    # are done; otherwise, we add it to the stack.
                    output = targets.get( adjacent, None )
                    if output is not None:
                        return output
                    else:
                        stack.append(adjacent)

        # End of loop. Move on to the next embedding of edge e.

    # If the search terminates without finding the ideal edge, then the ideal
    # edge must belong to a component that gets destroyed.
    return None


def printIdealEdges(surfaces):
    """
    Prints details of all ideal edges obtained by crushing annuli in the
    given list of normal surfaces.
    """
    for i, surf in enumerate(surfaces):
        if isAnnulus(surf):
            print( i, idealEdge(surf) )


def crushAnnuli(surfaces):
    """
    Crushes all annuli in the given packet of normal surfaces, and adds a
    Container of the resulting triangulations as a child of the given packet.
    """
    results = Container( "Crush annuli" )
    surfaces.insertChildLast(results)
    annulusCount = 0
    for surfNum, surf in enumerate(surfaces):
        if not isAnnulus(surf):
            continue
        annulusCount += 1

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
        idEdge = idealEdge(surf)
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
            if not comp.isValid():
                comp.setLabel( comp.label() + ": INVALID" )

                # Fill in invalid boundary.
                filled = PacketOfTriangulation3(comp)
                filled.finiteToIdeal()
                filled.intelligentSimplify()
                filledLabel = comp.adornedLabel(
                        "Closed by \"making ideal\"" )
                comp.insertChildLast(filled)

                # Decompose the filled manifold into prime pieces.
                summands = filled.summands()
                if len(summands) == 0:
                    filled.setLabel( filledLabel + ": S3" )
                elif len(summands) == 1:
                    # Try to combinatorially recognise this filled manifold.
                    std = StandardTriangulation.recognise( summands[0] )
                    if std is None:
                        name = "Not recognised"
                    else:
                        name = std.manifold().name()
                    filled.setLabel( filledLabel + ": {}".format(name) )
                else:
                    filled.setLabel( filledLabel + ": Non-prime" )
                    for sumNum, s in enumerate(summands):
                        summand = PacketOfTriangulation3(s)
                        filled.insertChildLast(summand)

                        # Try to combinatorially recognise this summand.
                        std = StandardTriangulation.recognise(summand)
                        if std is None:
                            name = "Not recognised"
                        else:
                            name = std.manifold().name()
                        summand.setLabel( "Summand #{}: {}".format(
                            sumNum, name ) )
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

                    # Try to recognise the drilled manifold.
                    if drilled.isSolidTorus():
                        name = "Ideal solid torus"
                    else:
                        # Try to combinatorially recognise after truncating
                        # the ideal vertex.
                        trunc = PacketOfTriangulation3(drilled)
                        drilled.insertChildLast(trunc)
                        trunc.idealToFinite()
                        trunc.intelligentSimplify()
                        std = StandardTriangulation.recognise(trunc)
                        if std is None:
                            name = "Not recognised"
                        else:
                            name = std.manifold().name()
                        trunc.setLabel( drilled.adornedLabel(
                            "Truncated" ) + ": {}".format(name) )
                    drilled.setLabel(
                            drilled.label() + ": {}".format(name) )

                # Decompose this component into prime pieces.
                summands = comp.summands()
                if len(summands) == 0:
                    comp.setLabel( comp.label() + ": S3" )
                elif len(summands) == 1:
                    # Try to combinatorially recognise this component.
                    std = StandardTriangulation.recognise( summands[0] )
                    if std is None:
                        name = "Not recognised"
                    else:
                        name = std.manifold().name()
                    comp.setLabel( comp.label() + ": {}".format(name) )
                else:
                    comp.setLabel( comp.label() + ": Non-prime" )
                    if compNum == idComp:
                        sumContainer = Container(
                                "Summands (Total {})".format(
                                    len(summands) ) )
                        comp.insertChildLast(sumContainer)
                    else:
                        sumContainer = comp
                    for sumNum, s in enumerate(summands):
                        summand = PacketOfTriangulation3(s)
                        sumContainer.insertChildLast(summand)

                        # Try to combinatorially recognise this summand.
                        std = StandardTriangulation.recognise(summand)
                        if std is None:
                            name = "Not recognised"
                        else:
                            name = std.manifold().name()
                        summand.setLabel( "Summand #{}: {}".format(
                            sumNum, name ) )
    results.setLabel( results.adornedLabel(
        "Total {}".format(annulusCount) ) )
