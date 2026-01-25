"""
Auxiliary functions for calculating how tetrahedra would get renumbered after
some tetrahedra get deleted.
"""


def tetRenumbering(doomed):
    """
    Returns a list describing how tetrahedra get renumbered after deleting the
    given doomed tetrahedra.

    In detail, letting r denote the returned list, for each index i of a
    tetrahedron T that is not doomed, r[i] will be the new index of T after
    all the doomed tetrahedra have been deleted. For an index i of a doomed
    tetrahedron, r[i] will be None.

    Precondition:
    --> doomed is nonempty, and contains only tetrahedra that all belong to
        the same 3-manifold triangulation.
    """
    tri = doomed[0].triangulation()
    doomedIndices = { tet.index() for tet in doomed }
    return tetIndexRenumbering( tri, doomedIndices )


def tetIndexRenumbering( tri, doomedIndices ):
    """
    Returns a list describing how tetrahedron in the given triangulation get
    renumbered after deleting the tetrahedra with the given doomed indices.

    In detail, letting r denote the returned list, for each index i of a
    tetrahedron T that is not doomed, r[i] will be the new index of T after
    all the doomed tetrahedra have been deleted. For a doomed index, i, r[i]
    will be None.
    """
    renum = []
    tetShift = 0
    for tetIndex in range( tri.size() ):
        if tetIndex in doomedIndices:
            renum.append(None)
            tetShift += 1
        else:
            renum.append( tetIndex - tetShift )
    return renum
