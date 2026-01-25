"""
Operations on triangulations that involve inserting a tetrahedron.
"""
from regina import *


def snapEdge( edge, check=True, perform=True ):
    """
    If the endpoints of the given edge are distinct and not both boundary,
    then uses a snapped ball to pinch these two endpoints together.

    This operation is equivalent to performing the following two operations:
    (1) Pinching the edge, which introduces a two-tetrahedron gadget with a
        single degree-one edge e at its heart.
    (2) Performing a 2-1 edge move on e.

    If check is True (the default), then this routine will check whether
    snapping the given edge is legal; otherwise, this routine will proceed
    under the assumption that the move is already known to be legal. If
    perform is True (the default), then this routine will actually perform
    the snap edge move if it has determined or assumed that the move is
    legal; otherwise, the triangulation containing the given edge will be
    left unchanged.

    If the triangulation containing the given edge is currently oriented,
    then this operation will preserve the orientation.

    Parameters:
    --> edge    The edge whose endpoints should be snapped together.

    Returns:
        True if and only if snapping the given edge is legal.
    """
    if check:
        # Endpoints need to be distinct and not both boundary.
        u = edge.vertex(0)
        v = edge.vertex(1)
        if u == v:
            return False
        if u.isBoundary() and v.isBoundary():
            return False
    if not perform:
        return True

    # Start by pinching the given edge.
    tri = edge.triangulation()
    tri.pinchEdge(edge)

    # To find the degree-one edge at the heart of the pinch edge gadget, look
    # at the last two tetrahedra in tri.
    found = False
    for tetIndex in [ tri.size() - 1, tri.size() - 2 ]:
        for edgeNum in range(6):
            e = tri.tetrahedron(tetIndex).edge(edgeNum)
            if e.degree() == 1:
                found = True
                break
        if found:
            break

    # Finish up by performing a 2-1 move on e.
    if not tri.twoOneMove( e, 0 ):
        if not tri.twoOneMove( e, 1 ):
            # This should never happen.
            raise RuntimeError( "Snap edge failed unexpectedly." )
    return True


def layerOn(edge):
    """
    Performs a layering upon the given boundary edge of a 3-manifold
    triangulation.

    This is almost equivalent to edge.triangulation().layerOn(edge). The only
    difference is that unlike Regina's implementation, this routine
    guarantees to produce an oriented triangulation if the original
    edge.triangulation() is oriented.

    Like Regina's layerOn() routine, this routine returns the newly layered
    tetrahedron T, and the new boundary edge created by the layering will be
    edge 5 of T (i.e., the edge joining vertices 2 and 3 of T).

    Pre-condition:
    --> The given edge is a boundary edge of a 3-manifold triangulation.
    --> The boundary triangles on either side of the given edge are distinct.

    Exceptions:
    --> InvalidArgument     The pre-conditions above do not hold. That is,
                            either the given edge is non-boundary, or the
                            same boundary triangle lies on both sides of it.

    Parameters:
    --> edge    The boundary edge upon which to layer.

    Returns:
        The new tetrahedron provided by the layering.
    """
    if not edge.isBoundary():
        raise InvalidArgument( "layerOn() requires a boundary edge." )
    fTet = edge.front().tetrahedron()
    fVer = edge.front().vertices()
    bTet = edge.back().tetrahedron()
    bVer = edge.back().vertices()
    if fTet.triangle( fVer[3] ) == bTet.triangle( bVer[2] ):
        raise InvalidArgument( "layerOn() requires an edge between two " +
                              "distinct boundary triangles." )

    # Up to rotation and/or reflection, the vertices of the front boundary
    # triangle are labelled as shown below left, and the vertices of the back
    # boundary triangle are labelled as shown below right.
    #
    #                     fVer[0]• •bVer[0]
    #                           /| |\
    #                          / | | \
    #                  fVer[2]•  | |  •bVer[3]
    #                          \ | | /
    #                           \| |/
    #                     fVer[1]• •bVer[1]
    #
    # If the triangulation is oriented, then the permutations fVer and bVer
    # will both have sign +1, so we can ensure that the triangulation remains
    # oriented after layering if we attach the new tetrahedron as follows:
    #
    #                             0
    #                             •
    #                            /|\
    #                           / | \
    #                         3•-----•2
    #                           \ | /
    #                            \|/
    #                             •
    #                             1
    #
    newTet = edge.triangulation().newTetrahedron()
    newTet.join( 2, fTet, fVer * Perm4(2,3) )
    newTet.join( 3, bTet, bVer * Perm4(2,3) )
    return newTet
