"""
Find non-surviving triangular orbits by traversing wedge cells.
"""
from enum import Enum, auto
from regina import *
from aux.quad import tetQuadType


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
        # at the *front* in the figure below, then the vertex numbers of tet
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
        eOrder = [ tet.edgeMapping(quadType), tet.edgeMapping(5-quadType) ]
        wedgeAdjacencies[teti] = dict()
        for i in range(2):
            # Encode wedge cell i via the following triple:
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

    # Count non-surviving triangular orbit types by explicitly traversing all
    # wedge cells.
    orbitCounts = { OrbitType.BOUNDARY: 0,
                   OrbitType.TRIVIAL_CYCLE: 0,
                   OrbitType.TWIST_PLUS: 0,
                   OrbitType.TWIST_MINUS: 0 }
    while wedgePerms:
        startWedge, endPerm = wedgePerms.popitem()
        startTeti, forwardFace, backwardFace = startWedge

        # Traverse *forwards* until one of the following occurs:
        #   --> We return to the start (in which case we have found a cycle
        #       of wedge cells).
        #   --> We hit the boundary.
        #   --> We reach a tetrahedron with no wedge cells.
        currentTet = tri.tetrahedron(startTeti)
        currentFace = forwardFace
        vertPerm = Perm4()
        forwardsEndedAtBoundary = False # Until proven otherwise.
        foundCycle = False
        while True:
            # Traverse across face gluing.
            adjTet = currentTet.adjacentTetrahedron(currentFace)
            if adjTet is None:
                forwardsEndedAtBoundary = True
                break
            adjTeti = adjTet.index()
            adjFace = currentTet.adjacentFace(currentFace)
            adjGluing = currentTet.adjacentGluing(currentFace)
            vertPerm = adjGluing * vertPerm

            # Have we reached another wedge cell?
            if adjTeti not in wedgeAdjacencies:
                # Forwards traversal ended in the interior, at a tetrahedron
                # with no wedge cells.
                break
            currentFace, adjWedge = wedgeAdjacencies[adjTet.index()][adjFace]

            # Have we returned to the start? If so, then use the twistSign
            # dictionary to determine which direction this cycle twists.
            if adjWedge == startWedge:
                vertPerm = endPerm * vertPerm
                cycleSign = twistSign[startWedge][ vertPerm[backwardFace] ]
                if cycleSign == 1:
                    orbitCounts[ OrbitType.TWIST_PLUS ] += 1
                elif cycleSign == -1:
                    orbitCounts[ OrbitType.TWIST_MINUS ] += 1
                else:   # cycleSign == 0
                    orbitCounts[ OrbitType.TRIVIAL_MINUS ] += 1
                foundCycle = True
                break

            # Traverse across the new wedge cell.
            # No need to set currentFace since that was done earlier.
            vertPerm = wedgePerms.pop(adjWedge) * vertPerm
            currentTet = adjTet
        # End of forwards traversal.
        if foundCycle:
            continue

        # At this point, we definitely do not have a cycle of wedge cells,
        # but we still need to traverse backwards to both:
        #   (a) find out whether we have a (non-surviving) triangular orbit
        #       where both ends are at the boundary; and
        #   (b) remove all wedges in this orbit from wedgePerms.
        currentTet = tri.tetrahedron(startTeti)
        currentFace = backwardFace
        while True:
            # Traverse across face gluing.
            adjTet = currentTet.adjacentTetrahedron(currentFace)
            if adjTet is None:
                # Backwards traversal has ended at boundary.
                if forwardsEndedAtBoundary:
                    orbitCounts[ OrbitType.BOUNDARY ] += 1
                break
            adjTeti = adjTet.index()
            adjFace = currentTet.adjacentFace(currentFace)
            adjGluing = currentTet.adjacentGluing(currentFace)

            # Have we reached another wedge cell?
            if adjTeti not in wedgeAdjacencies:
                # Backwards traversal ended in the interior, at a tetrahedron
                # with no wedge cells.
                break
            currentFace, adjWedge = wedgeAdjacencies[adjTet.index()][adjFace]
            currentTet = adjTet

            # Don't forget to remove from wedgePerms, so that we only
            # traverse each wedge exactly once.
            wedgePerms.pop(adjWedge)
        # End of backwards traversal.

    # All done!
    return orbitCounts
