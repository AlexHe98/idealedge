"""
Traverse wedge cells to detect lost L(3,1) components.

In the context of a Seifert fibre space with (some component of the) boundary
given by an ideal loop, such lost components correspond to a lost fibre of
multiplicity 3.
"""
from regina import *
from quadtype import tetQuadType


def wedgeLoops(surf):
    """
    Detects loops of wedge cells induced by the given normal surface.

    This routine returns a list of such loops, each of which is encoded as a
    pair (R,T), where:
    --> R is a single wedge cell in the loop, chosen as a representative for
        the entire loop.
    --> T is 0 if the loop has no twist, and either +1 or -1 if the loop has
        a twist.

    We have exactly two wedge cells per tetrahedron intersecting surf in a
    positive number of quads. Each wedge cell is encoded as a pair (i,s),
    where:
    --> i is the index of the tetrahedron containing the wedge cell; and
    --> s is 0 if the wedge cell is incident to edge q of tetrahedron i,
        where q is the quad type, and 1 if the wedge cell is incident to edge
        5-q of tetrahedron i.

    The sign of the twist is determined relative to the vertex labelling of
    the tetrahedron containing the representative wedge cell. Thus, in an
    oriented triangulation, wedge loops with the same sign will twist in the
    same direction.
    """
    tri = surf.triangulation()

    # Find all wedge cells.
    wedgePerms = dict()         # Map wedge cells to vertex permutations.
    wedgeAdjacencies = dict()
    twistSign = dict()
    for tet in tri.tetrahedra():
        teti = tet.index()
        quadType = tetQuadType( surf, teti )
        if quadType is None:
            continue

        # We have at least one quad in this tet, and hence two wedge cells.
        #
        # If we index the wedge cells by i in {0,1}, and if wedge i is drawn
        # at the front in the figure below, then the vertex numbers of tet
        # will be as shown in the figure.
        #
        #               eOrder[1-i][1]
        #                      •
        #                     /|\
        #                    / | \
        #                   /__|__\
        #                  /|  |  |\
        #                 / |  ↑  | \
        #   eOrder[i][0] •--|-→|→-|--• eOrder[i][1]
        #                 \ |  ↑  | /
        #                  \|__|__|/
        #                   \  |  /
        #                    \ | /
        #                     \|/
        #                      •
        #               eOrder[1-i][0]
        #
        eOrder = [ tet.edgeMapping(quadType), tet.edgeMapping(5-quadType) ]
        wedgeAdjacencies[teti] = dict()
        for i in range(2):
            wedge = ( teti, eOrder[i][0], eOrder[i][1] )
            wedgePerms[wedge] = eOrder[i] * Perm4(0,1) * eOrder[i].inverse()

            # Later on, we might need to compute the sign of a twist at this
            # wedge cell. For this, we will use eOrder[i][1] as a reference
            # vertex, and determine the sign by examining where this
            # reference vertex gets mapped to.
            twistSign[wedge] = {
                    eOrder[i][1]: 0,
                    eOrder[i][3]: 1,
                    eOrder[i][2]: -1 }

            # Record adjacency of faces across wedge cells; we need this
            # information to traverse across wedge cells.
            for ii in range(2):
                wedgeAdjacencies[teti][eOrder[i][ii]] =\
                        ( eOrder[i][1-ii], wedge )

    # Traverse all wedge cells.
    loops = set()
    while wedgePerms:
        startWedge, endPerm = wedgePerms.popitem()
        startTeti, currentFace, referenceVertex = startWedge
        currentTet = tri.tetrahedron(startTeti)
        vertPerm = Perm4()

        # Traverse until one of the following occurs:
        #   --> We return to the start (in which case we have found a loop of
        #       wedge cells).
        #   --> We reach a wedge cell that we already previously traversed
        #       (in which case we do not have a loop of wedge cells).
        #   --> We reach a tetrahedron with no wedge cells (in which case we
        #       again do not have a loop of wedge cells).
        while True:
            # Traverse across face gluing.
            adjTet = currentTet.adjacentTetrahedron(currentFace)
            if adjTet is None:
                # No adjacent tet, so definitely not traversing a loop of
                # wedge cells.
                break
            adjTeti = adjTet.index()
            adjFace = currentTet.adjacentFace(currentFace)
            adjGluing = currentTet.adjacentGluing(currentFace)
            vertPerm = adjGluing * vertPerm

            # Have we reached a new wedge cell?
            if adjTeti not in wedgeAdjacencies:
                # No wedge cells in adjTet, and hence we are not traversing
                # a loop of wedge cells.
                break
            currentFace, adjWedge = wedgeAdjacencies[adjTet.index()][adjFace]

            # Have we already previously traversed the new wedge cell? If so,
            # then either:
            #   --> we have returned to the start of a loop of wedge cells;
            #       or
            #   --> we are not traversing a loop of wedge cells.
            if adjWedge not in wedgePerms:
                # If we are back to the start of a loop, then use the
                # twistSign dictionary to determine which direction this loop
                # twists.
                if adjWedge == startWedge:
                    vertPerm = endPerm * vertPerm
                    loopSign = twistSign[startWedge][
                            vertPerm[referenceVertex] ]
                    loops.add( ( startWedge, loopSign ) )

                # Regardless of whether or not we had a loop, there is no
                # further traversal we can do.
                break

            # Traverse across the new wedge cell.
            # No need to set currentFace since that was done earlier.
            vertPerm = wedgePerms.pop(adjWedge) * vertPerm
            currentTet = adjTet

    # All done!
    return loops
