"""
Pinch an ideal loop without forgetting the meridian.
"""
from regina import *


def pinch(loop):
    """
    Returns the meridian of the given ideal loop, represented as a curve in a
    vertex link of a newly constructed ideal triangulation given by pinching
    the given ideal loop.

    This routine returns a triple (t, v, f), where:
    --> t is a tetrahedron incident to the meridian;
    --> v is a vertex number of t whose corresponding vertex-linking triangle
        is incident to the meridian; and
    --> f is a triangle number of t whose corresponding triangle is incident
        to both the meridian and vertex v.

    If the triangulation containing the given ideal loop is oriented, then
    the same is guaranteed to be true for the new ideal triangulation.

    Pre-condition:
    --> The given ideal loop has length one.
    """
    if len(loop) != 1:
        raise ValueError("Can only pinch an ideal loop of length one.")
    tri = Triangulation3( loop.triangulation() )

    # Before we make any modifications to tri, we should find the face along
    # which we will insert the pinch gadget. Specifically, the face will be
    # tet.triangle( ver[3] ).
    emb = tri.edge( loop[0] ).embedding(0)
    tet = emb.tetrahedron()
    ver = emb.vertices()

    # Introduce a pinched 3-ball such that the core curve of the annulus link
    # can be found at the triangle linking vertex merVert of pinchTet[0],
    # along the edge running through face merFace of pinchTet[0].
    pinchTet, merVert, merFace = _addPinchGadgetWithMeridian(tri)

    # Record the face gluing before undoing it.
    glue = tet.adjacentGluing( ver[3] )
    adj = tet.adjacentTetrahedron( ver[3] )
    tet.unjoin( ver[3] )

    # Now insert the pinched 3-ball.
    #
    # By construction, the pinched 3-ball has boundary faces given by faces 0
    # and 2 of pinchedTet[1], such that:
    #   --> in face 0, vertices 2 and 3 are pinched together
    #   --> in face 2, vertices 0 and 3 are pinched together
    #
    #   Vertices and edges of pinchTet[1]
    #        0*---*1
    #         |  /
    #         | /
    #         |/
    #      2~3*----loop wrapped around equator of pinchTet[0]
    #   (vertex 2 "in front")
    #
    # Note that if the original triangulation was oriented, then the ver
    # permutation is guaranteed to have sign +1. This ensures that all the
    # new gluings have sign -1, and hence that the orientation is preserved
    # after inserting the pinched 3-ball.
    pinchTet[1].join( 0, tet, ver * Perm4(3,2,0,1) )
    pinchTet[1].join( 2, adj, glue * ver * Perm4(0,2,3,1) )

    # All done!
    return pinchTet[0], merVer, merFace


def _addPinchGadgetWithMeridian(tri):
    pinchTet = tri.newTetrahedra(2)

    # Form a snapped 3-ball with pinchTet[0] so that:
    #   --> the northern hemisphere is given by triangle 0
    #   --> the southern hemisphere is given by triangle 1
    #
    #   Vertices and edges of pinchTet[0]
    #        1*
    #         |
    #         |
    #         |
    #      2~3*----equator loop
    #         |
    #         |
    #         |
    #        0*
    #   (vertex 2 "in front")
    #
    pinchTet[0].join( 2, pinchTet[0], Perm4(2,3) )

    # Attach pinchTet[1] to the northern hemisphere.
    #
    #   Vertices and edges of pinchTet[1]
    #        0*---*1
    #         |  /
    #         | /
    #         |/
    #      2~3*----loop wrapped around equator of pinchTet[0]
    #   (vertex 2 "in front")
    #
    pinchTet[0].join( 0, pinchTet[1], Perm4(0,1) )

    # To turn this into a pinched 3-ball, the idea is to wrap face 012 of
    # pinchTet[1] around the southern hemisphere, such that edge 02 goes
    # around the equator.
    pinchTet[0].join( 1, pinchTet[1], Perm4(1,3,0,2) )

    # The result is a pinched 3-ball with boundary faces given by faces 0 and
    # 2 of pinchedTet[1], such that:
    #   --> in face 0, vertices 2 and 3 are pinched together
    #   --> in face 2, vertices 0 and 3 are pinched together
    #
    #   loop0~2   *1
    #         |  /
    #         | /
    #         |/
    #    0~2~3*----loop2~3
    #
    # If we cut open a face F of the given triangulation and glue in this
    # pinched 3-ball, the result will be to pinch (a curve parallel to) an
    # edge of F down to a point.
    #
    # The meridian of the pinched curve appears as a single edge in the link
    # of the pinched vertex. Specifically, if we consider the triangle L
    # linking vertex 3 of pinchTet[0], the meridian is given by the edge of L
    # that runs through face 0 of pinchTet[0].
    return pinchTet, 3, 0


if __name__ == "__main__":
    tri = Triangulation3()
    _addPinchGadgetWithMeridian(tri)
    print( tri.detail() )
    print()
    for v in tri.vertices():
        if not v.isValid():
            print( v.index() )
