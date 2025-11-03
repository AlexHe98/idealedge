"""
Helper routines for accessing quads and quad types of a normal surface.
"""


def tetQuadType( surf, tetIndex ):
    """
    Returns the quad type in which the given normal surface intersects the
    tetrahedron with the given index, or None if there is no such quad.
    """
    quads = tetQuads( surf, tetIndex )
    if quads is None:
        return None
    else:
        return quads[0]


def tetQuads( surf, tetIndex ):
    """
    Returns the quad type and the number of quads in which the given normal
    surface intersects the tetrahedron with the given index, or None if there
    is no such quad.
    """
    for quadType in range(3):
        quadCount = surf.quads( tetIndex, quadType ).safeLongValue()
        if quadCount > 0:
            return ( quadType, quadCount )
    return None
