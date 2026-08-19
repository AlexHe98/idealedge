"""
Find non-surviving triangular orbits by traversing wedge cells.
"""
from enum import Enum, auto
from regina import *
from aux.quad import tetQuadType
#TODO What about non-surviving triangular orbits which go from boundary to
#   boundary (and hence give deleted 3-ball components)?


class NonSurvivingTriangularOrbitType(Enum):
    """
    An enumeration of types of non-surviving triangular orbits.

    In detail, we have the following possible types:
    --> BOUNDARY        A triangular orbit whose core fibre forms an arc with
                        both endpoints lying in real boundary.
    --> TRIVIAL_CYCLE   A triangular orbit whose core fibre forms a loop with
                        no twist.
    --> TWIST_PLUS      A triangular orbit whose core fibre forms a loop with
                        a twist. The sign of the twist is +1, as determined
                        by the orientation of some tetrahedron incident to
                        the orbit.
    --> TWIST_MINUS     A triangular orbit whose core fibre forms a loop with
                        a twist. The sign of the twist is -1, as determined
                        by the orientation of some tetrahedron incident to
                        the orbit.
    As defined, the signs of the twists depend on the chosen tetrahedron.
    Thus, these signs are most meaningful when we are working with an
    oriented triangulation.
    """
    BOUNDARY = auto()
    TRIVIAL_CYCLE = auto()
    TWIST_PLUS = auto()
    TWIST_MINUS = auto()
    pass


def nonSurvivingTriangularOrbitCounts(surf):
    """
    Counts the number of non-surviving triangular orbits of each possible
    type.

    In detail, this routine returns a dictionary N with keys given by the
    NonSurvivingTriangularOrbitType enumeration, such that for each key k,
    N[k] counts the number of non-surviving triangular orbits of type k.

    Pre-condition:
    --> surf.triangulation() is oriented.
    """
    OrbitType = NonSurvivingTriangularOrbitType
    tri = surf.triangulation()

    # Find all wedge cells.
    #
    # We have exactly two wedge cells per tetrahedron intersecting surf in a
    # positive number of quads. Each wedge cell is encoded as a pair (i,s),
    # where:
    #   --> i is the index of the tetrahedron containing the wedge cell; and
    #   --> s is 0 if the wedge cell is incident to edge q of tetrahedron i,
    #       where q is the quad type, and 1 if the wedge cell is incident to
    #       edge 5-q of tetrahedron i.
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

    #TODO Need to track wedge cells incident to boundary.

    # Traverse all wedge cells.
    orbitCounts = { OrbitType.BOUNDARY: 0,
                   OrbitType.TRIVIAL_CYCLE: 0,
                   OrbitType.TWIST_PLUS: 0,
                   OrbitType.TWIST_MINUS: 0 }
    while wedgePerms:
        startWedge, endPerm = wedgePerms.popitem()
        startTeti, currentFace, referenceVertex = startWedge
        currentTet = tri.tetrahedron(startTeti)
        vertPerm = Perm4()

        # Traverse until one of the following occurs:
        #   --> We return to the start (in which case we have found a cycle
        #       of wedge cells).
        #   --> We reach a wedge cell that we already previously traversed
        #       (in which case we do not have a cycle of wedge cells).
        #   --> We reach a tetrahedron with no wedge cells (in which case we
        #       again do not have a cycle of wedge cells).
        while True:
            # Traverse across face gluing.
            adjTet = currentTet.adjacentTetrahedron(currentFace)
            if adjTet is None:
                # No adjacent tet, so definitely not traversing a cycle of
                # wedge cells.
                break
            adjTeti = adjTet.index()
            adjFace = currentTet.adjacentFace(currentFace)
            adjGluing = currentTet.adjacentGluing(currentFace)
            vertPerm = adjGluing * vertPerm

            # Have we reached a new wedge cell?
            if adjTeti not in wedgeAdjacencies:
                # No wedge cells in adjTet, and hence we are not traversing
                # a cycle of wedge cells.
                break
            currentFace, adjWedge = wedgeAdjacencies[adjTet.index()][adjFace]

            # Have we already previously traversed the new wedge cell? If so,
            # then either:
            #   --> we have returned to the start of a cycle of wedge cells;
            #       or
            #   --> we are not traversing a cycle of wedge cells.
            if adjWedge not in wedgePerms:
                # If we are back to the start of a cycle, then use the
                # twistSign dictionary to determine which direction this
                # cycle twists.
                if adjWedge == startWedge:
                    vertPerm = endPerm * vertPerm
                    cycleSign = twistSign[startWedge][
                            vertPerm[referenceVertex] ]
                    if cycleSign == 1:
                        orbitCounts[ OrbitType.TWIST_PLUS ] += 1
                    elif cycleSign == -1:
                        orbitCounts[ OrbitType.TWIST_MINUS ] += 1
                    else:   # cycleSign == 0
                        orbitCounts[ OrbitType.TRIVIAL_MINUS ] += 1

                # Regardless of whether or not we had a cycle, there is no
                # further traversal we can do.
                break

            # Traverse across the new wedge cell.
            # No need to set currentFace since that was done earlier.
            vertPerm = wedgePerms.pop(adjWedge) * vertPerm
            currentTet = adjTet

    # All done!
    return orbitCounts
