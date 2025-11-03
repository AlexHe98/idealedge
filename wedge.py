"""
Traverse wedge cells to detect lost L(3,1) components.
"""
from sys import argv
from regina import *
from parallel import tetQuadType


def wedgeLoops(surf):
    """
    Detects loops of wedge cells induced by the given normal surface.
    """
    tri = surf.triangulation()

    # Find all wedge cells.
    wedgePerms = dict()         # Map wedge cells to vertex permutations.
    wedgeAdjacencies = dict()
    for teti in range( tri.size() ):
        quadType = tetQuadType( surf, teti )
        if quadType is None:
            continue

        # We have at least one quad in this tet, and hence two wedge cells.
        #
        # If we index the wedge cells by i in {0,1}, and if wedge i is drawn
        # at the front in the figure below, then the vertex numbers of tet
        # will be as shown in the figure.
        #
        #               eOrder[1-i][0]
        #                      *
        #                     /|\
        #                    / | \
        #                   /__|__\
        #                  /|  |  |\
        #                 / |  |  | \
        #   eOrder[i][0] *--|--|--|--* eOrder[i][1]
        #                 \ |  |  | /
        #                  \|__|__|/
        #                   \  |  /
        #                    \ | /
        #                     \|/
        #                      *
        #               eOrder[1-i][1]
        #
        eOrder = [ Edge3.ordering(quadType), Edge3.ordering(5-quadType) ]
        wedgeAdjacencies[teti] = dict()
        for i in range(2):
            wedge = ( teti, eOrder[i][0], eOrder[i][1] )
            wedgePerms[wedge] = eOrder[i] * Perm4(0,1) * eOrder[i].inverse()
            for ii in range(2):
                # Record adjacency of faces across wedge cells; we need this
                # information to traverse across wedge cells.
                wedgeAdjacencies[teti][eOrder[i][ii]] =\
                        ( eOrder[i][1-ii], wedge )

    # Traverse all wedge cells.
    #TODO Decide on returned data structure.
    #   Currently using a set of pairs of the form (r,t), where r is a
    #   representative wedge cell, and t is True if and only if the loop has
    #   a twist.
    loops = set()
    while wedgePerms:
        startWedge, endPerm = wedgePerms.popitem()
        startTeti, currentFace, endFace = startWedge
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
                # If we are back to the start of a loop, then make a note of
                # whether or not this loop has a twist.
                if adjWedge == startWedge:
                    vertPerm = endPerm * vertPerm
                    loops.add( ( startWedge, not vertPerm.isIdentity() ) )

                # Regardless of whether or not we had a loop, there is no
                # further traversal we can do.
                break

            # Traverse across the new wedge cell.
            # No need to set currentFace since that was done earlier.
            vertPerm = wedgePerms[adjWedge] * vertPerm
            currentTet = adjTet

    # All done!
    return loops


if __name__ == "__main__":
    sig = argv[1]
    num = int( argv[2] )
    tri = Triangulation3.fromIsoSig(sig)
    qvsurfs = NormalSurfaces( tri, NormalCoords.NS_QUAD )
    surf = qvsurfs[num]

    #TODO TEST
    print( wedgeLoops(surf) )
