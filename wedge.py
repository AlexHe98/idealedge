"""
Traverse wedge cells to detect lost L(3,1) components.
"""
from regina import *
from parallel import tetQuadType


def wedgeLoops(surf):
    """
    Detects loops of wedge cells induced by the given normal surface.
    """
    #TODO Decide on returned data structure.
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
        #                  eOrder[i][0]
        #                        *
        #                       /|\
        #                      / | \
        #                     /__|__\
        #                    /|  |  |\
        #                   / |  |  | \
        #   eOrder[1-i][0] *--|--|--|--* eOrder[1-i][1]
        #                   \ |  |  | /
        #                    \|__|__|/
        #                     \  |  /
        #                      \ | /
        #                       \|/
        #                        *
        #                  eOrder[i][1]
        #
        eOrder = [ Edge3.ordering(quadType), Edge3.ordering(5-quadType) ]
        wedgeAdjacencies[teti] = dict()
        for i in range(2):
            wedge = ( teti, eOrder[i][0], eOrder[i][1] )
            wedgePerms[wedge] = eOrder[i] * Perm4(2,3) * eOrder[i].inverse()
            for ii in range(2):
                # Record adjacency of faces across wedge cells; we need this
                # information to traverse across wedge cells.
                wedgeAdjacencies[teti][eOrder[1-i][ii]] =\
                        ( eOrder[1-i][1-ii], wedge )

    # Traverse all wedge cells.
    while wedgePerms:
        startWedge, startPerm = wedgePerms.popitem()
        startTeti, currentFace, endFace = startWedge
        currentTet = tri.tetrahedron(startTeti)
        currentPerm = Perm4()

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
            currentPerm = adjGluing * currentPerm

            # Have we reached a new wedge cell?
            if adjTeti not in wedgeAdjacencies:
                # No wedge cells in adjTet.
                #TODO
                break
            nextFace, adjWedge = wedgeAdjacencies[adjTet.index()][adjFace]
            #TODO

            # Traverse across adjacent wedge cell.
            #TODO
            pass
    #TODO
    raise NotImplementedError()
