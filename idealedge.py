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

    Specifically, this routine returns a pair (t,e), where t is the index
    (after crushing) of a tetrahedron that will meet the ideal edge after
    crushing, and e is an edge number of this tetrahedron that corresponds to
    the ideal edge.

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
        wt = annulus.edgeWeight(ei)
        for s in range( wt + 1 ):
            segments.append( (ei,s) )

    # Find "target" segments.
    targets = dict()
    for tet in tri.tetrahedra():
        teti = tet.index()
        hasQuads = False
        for q in range(3):
            if annulus.quads( teti, q ) > 0:
                hasQuads = True
                break
        if hasQuads:
            continue

        # No quads in tet, so there is a cell in the centre that survives
        # crushing. Find the edges of this cell that survive.
        for en in range(6):
            v = tet.edgeMapping(en)[0]
            ei = tet.edge(en).index()
            s = annulus.triangles( teti, v )
            targets[ (ei,s) ] = ( teti, en )

    # Starting segment for depth-first search.
    for e in tri.edges():
        ei = e.index()
        if ( not e.isBoundary() ) or ( annulus.edgeWeight(ei) <= 1 ):
            continue
        start = ( ei, 1 )

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
        wt = annulus.edgeWeight(ei)
        visited.add(current)    # Record as visited now, so we don't forget.
        for emb in e.embeddings():
            tet = emb.tetrahedron()
            teti = tet.index()
            en = emb.face()
            ver = emb.vertices()

            # To locate the relevant parallel cells in tet, need to get the
            # normal coordinates incident to e.
            f = [ annulus.triangles( teti, ver[i] ) for i in range(2) ]
            q = 0
            qType = None
            for qt in range(3):
                if qt in { en, 5-en }:
                    continue
                quads = annulus.quads( teti, qt )
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
                    verOther = tet.faceMapping(enOther)

                    # Is the current segment adjacent to a segment of the
                    # edge numbered enOther?
                    if verOther[0] == ver[0]:
                        adjacent = ( eiOther, seg )
                    elif verOther[1] == ver[0]:
                        wtOther = annulus.edgeWeight(eiOther)
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
                if ver[0] not in s0:
                    side0 = {0,1,2,3} - side0
                qDepth = seg - f[0]     # 1 <= qDepth <= q - 1
                for enOther in range(6):
                    if enOther in { en, qType, 5 - qType }:
                        continue
                    eiOther = tet.edge(enOther).index()
                    verOther = tet.faceMapping(enOther)

                    # The current segment is adjacent to a segment of the
                    # edge numbered enOther. Find this adjacent segment.
                    fOther = annulus.triangles( teti, verOther[0] )
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
                    verOther = tet.faceMapping(enOther)

                    # Is the current segment adjacent to a segment of the
                    # edge numbered enOther?
                    if verOther[0] == ver[1]:
                        adjacent = ( eiOther, wt - seg )
                    elif verOther[1] == ver[1]:
                        wtOther = annulus.edgeWeight(eiOther)
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

    # If the search terminates without finding the ideal edge, then something
    # must have gone wrong.
    raise RuntimeError( "Failed to find ideal edge." )
