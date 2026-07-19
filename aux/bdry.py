"""
Auxiliary routines for querying the boundary of 3-manifold triangulations.
"""
from regina import *


def hasOnlyMinimalRealTorusBoundaryComponents(tri):
    """
    Returns True if and only if every boundary component of tri is a real
    two-triangle torus.

    In particular, this routine vacuously returns True if tri is closed.
    """
    for bc in tri.boundaryComponents():
        if ( not ( bc.isReal() and bc.size() == 2 and bc.isOrientable() and
                  bc.eulerChar() == 0 ) ):
            return False
    return True
