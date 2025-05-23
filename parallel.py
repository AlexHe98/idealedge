"""
Find the ideal edges after crushing a normal surface.
"""
from sys import argv
from regina import *
from loop import NotLoop, IdealLoop


def findIdealEdges( surf, start, targets=None ):
    """
    Returns details of the ideal edge that corresponds to the given start
    segment after crushing surf.

    Specifically, if the ideal edge belongs to a component that gets
    destroyed after crushing, then this routine returns None. Otherwise, this
    routine returns a pair (t, e), where t is the index after crushing of a
    tetrahedron that will be incident to the ideal edge, and e is an edge
    number of this tetrahedron that corresponds to the ideal edge.

    If the dictionary of surviving segments has been precomputed using the
    _survivingSegments() routine, then this can be supplied using the
    optional targets argument. Otherwise, this routine will compute the
    surviving segments for itself.
    """
    #TODO Fix documentation/comments.
    tri = surf.triangulation()
    if targets is None:
        targets = _survivingSegments(surf)

    # If the start segment is one of the targets, then we are already done.
    found = set()
    output = targets.get( start, None )
    if output is not None:
        found.add(output)
        #return output

    # Otherwise, we find the ideal edge using depth-first search.
    stack = [start]
    visited = set()
    while stack:
        current = stack.pop()
        if current in visited:
            continue

        # We haven't visited the current segment yet, so we need to find all
        # segments that are adjacent to it along parallel cells or faces.
        ei, seg = current
        e = tri.edge(ei)
        wt = surf.edgeWeight(ei).safeLongValue()
        visited.add(current)    # Record as visited now, so we don't forget.
        for emb in e.embeddings():
            tet = emb.tetrahedron()
            teti = tet.index()
            en = emb.face()
            ver = emb.vertices()

            # To locate the relevant parallel cells and faces in tet, need to
            # get the normal coordinates incident to e.
            f = [ surf.triangles( teti, ver[i] ).safeLongValue()
                    for i in range(2) ]
            q = 0
            qType = None
            for qt in range(3):
                if qt in { en, 5-en }:
                    continue
                quads = surf.quads( teti, qt ).safeLongValue()
                if quads > 0:
                    q = quads
                    qType = qt
                    break

            # Does the current segment belong to a parallel cell or face in
            # this tet? If so, then we need to find all adjacent segments.
            if seg < f[0]:
                # The current segment belongs to a parallel triangular cell
                # at vertex ver[0].
                for otherEnd in range(4):
                    if otherEnd in { ver[0], ver[1] }:
                        continue

                    # The current segment is adjacent to a segment of the
                    # edge with endpoints ver[0] and otherEnd.
                    eiOther = tet.edge( ver[0], otherEnd ).index()
                    enOther = Edge3.edgeNumber[ver[0]][otherEnd]
                    verOther = tet.edgeMapping(enOther)
                    if verOther[0] == ver[0]:
                        adjacent = ( eiOther, seg )
                    else:
                        wtOther = surf.edgeWeight(eiOther).safeLongValue()
                        adjacent = ( eiOther, wtOther - seg )

                    # If the adjacent segment is one of the targets, then we
                    # are done; otherwise, we add it to the stack.
                    output = targets.get( adjacent, None )
                    if output is not None:
                        found.add(output)
                        #return output
                    else:
                        stack.append(adjacent)
            elif seg > f[0] + q:
                # The current segment belongs to a parallel triangular cell
                # at vertex ver[1].
                for otherEnd in range(4):
                    if otherEnd in { ver[0], ver[1] }:
                        continue

                    # The current segment is adjacent to a segment of the
                    # edge with endpoints ver[1] and otherEnd.
                    eiOther = tet.edge( ver[1], otherEnd ).index()
                    enOther = Edge3.edgeNumber[ver[1]][otherEnd]
                    verOther = tet.edgeMapping(enOther)
                    if verOther[0] == ver[1]:
                        adjacent = ( eiOther, wt - seg )
                    else:
                        wtOther = surf.edgeWeight(eiOther).safeLongValue()
                        adjacent = ( eiOther, wtOther - wt + seg )

                    # If the adjacent segment is one of the targets, then we
                    # are done; otherwise, we add it to the stack.
                    output = targets.get( adjacent, None )
                    if output is not None:
                        found.add(output)
                        #return output
                    else:
                        stack.append(adjacent)
            elif q > 0:
                # The quadrilaterals divide tet into two "sides". The edge
                # opposite this segment has endpoints lying on different
                # sides, so we label these endpoints accordingly.
                side = [ { 0, qType + 1 } ]
                side.append( {0,1,2,3} - side[0] )
                if ver[0] not in side[0]:
                    side[0], side[1] = side[1], side[0]
                side[0].remove( ver[0] )
                side[1].remove( ver[1] )
                opp = [ side[0].pop(), side[1].pop() ]

                # Find all edges containing segments that are adjacent to the
                # current segment.
                qDepth = seg - f[0]     # 0 <= qDepth <= q
                if qDepth == 0:
                    adjEndpoints = [ [ ver[0], opp[1] ] ]
                elif qDepth == q:
                    adjEndpoints = [ [ opp[0], ver[1] ] ]
                else:
                    adjEndpoints = [
                            [ ver[0], opp[1] ],
                            [ opp[0], ver[1] ],
                            [ opp[0], opp[1] ] ]
                for start, end in adjEndpoints:
                    # The current segment is adjacent to a segment of the
                    # edge going from start to end.
                    eiAdj = tet.edge( start, end ).index()
                    enAdj = Edge3.edgeNumber[start][end]
                    verAdj = tet.edgeMapping(enAdj)
                    triangles = surf.triangles(
                            teti, verAdj[0] ).safeLongValue()
                    if verAdj[0] == start:
                        adjacent = ( eiAdj, triangles + qDepth )
                    else:
                        adjacent = ( eiAdj, triangles + q - qDepth )

                    # If the adjacent segment is one of the targets, then we
                    # are done; otherwise, we add it to the stack.
                    output = targets.get( adjacent, None )
                    if output is not None:
                        found.add(output)
                        #return output
                    else:
                        stack.append(adjacent)

        # End of loop. Move on to the next embedding of edge e.

    # If the search terminates without finding the ideal edge, then the ideal
    # edge must belong to a component that gets destroyed.
    return found
    #return None


def _survivingSegments(surf):
    """
    Uses the given normal surface to divide the edges of the ambient
    triangulation into segments, and returns a dictionary describing the
    segments that would survive after crushing surf.

    In detail, the keys of the returned dictionary will be surviving segments
    encoded as pairs of the form (ei, s), where:
    --> ei is an edge index; and
    --> s is a segment number from 0 to w, inclusive, where w is the weight
        of surf on edge ei.
    The segments for each edge e are numbered in ascending order from the one
    incident to e.vertex(0) to the one incident to e.vertex(1). The returned
    dictionary will map each such segment to a pair (t, en), where:
    --> t is the index after crushing of a tetrahedron that will be incident
        to the segment in question; and
    --> en is an edge number (from 0 to 5, inclusive) of this tetrahedron
        that corresponds to the segment in question.
    """
    tri = surf.triangulation()
    survivors = dict()
    shift = 0
    for tet in tri.tetrahedra():
        teti = tet.index()
        hasQuads = False
        for q in range(3):
            if surf.quads( teti, q ).safeLongValue() > 0:
                hasQuads = True
                break
        if hasQuads:
            # Presence of quads means tet is destroyed by crushing, which
            # will shift all larger tetrahedron indices down by one.
            shift += 1
            continue

        # No quads in tet, so there is a cell in the centre that survives
        # crushing. Find the edges of this cell that survive.
        for en in range(6):
            v = tet.edgeMapping(en)[0]
            ei = tet.edge(en).index()
            s = surf.triangles( teti, v ).safeLongValue()
            survivors[ (ei,s) ] = ( teti - shift, en )

    # Done!
    return survivors


if __name__ == "__main__":
    ## SFS [D: (2,1) (2,1)] U/m SFS [D: (2,1) (2,1)], m = [ 0,-1 | 1,0 ] : #1
    #sig = "hLALMkbcceffggemkbtibj"
    sig = argv[1]
    num = int( argv[2] )
    tri = Triangulation3.fromIsoSig(sig)
    qvsurfs = NormalSurfaces( tri, NormalCoords.NS_QUAD )
    surf = qvsurfs[num]
    #print(surf)

    #TODO
    survivors = _survivingSegments(surf)
    found = []
    for e in tri.edges():
        ei = e.index()
        wt = surf.edgeWeight(ei).safeLongValue()

        # Find ideal edges.
        for s in range( 1, wt ):
            result = findIdealEdges( surf, ( ei, s ), survivors )
            if result not in found:
                found.append(result)
    print()
    for f in found:
        print(f)
