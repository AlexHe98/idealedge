"""
Randomisation for Regina's Triangulation3.

This bespoke implementation should eventually be replaced by
Triangulation3.simplifyUpDown(), which is expected to appear as part of the
release of Regina 8.0.
"""
from regina import *


def randomise(tri):
    """
    Attempts to randomly retriangulate the given triangulation.

    This routine works by performing lots of random 2-3 moves, before
    attempting to simplify the triangulation again. As such, this routine is
    often useful for escaping "wells" when tri.simplify() gets stuck.

    If tri is currently oriented, then this routine guarantees to preserve
    the orientation.

    Adapted from SnapPea's randomize_triangulation().
    """
    RandomEngine.reseedWithHardware()
    randomisation = 4       # Hard-coded value copied from SnapPea.
    origSize = tri.size()
    count = randomisation * origSize
    while count > 0:
        count -= 1

        # Attempt a random 2-3 move.
        if tri.pachner( tri.triangle( RandomEngine.rand(
            tri.countTriangles() ) ) ):
            # Try to force future random 2-3 moves to make "interesting"
            # changes.
            _simplifyBasic(tri)
            if tri.size() < origSize:
                # We already succeeded in escaping the well, so we might
                # as well terminate early.
                break

    # Finish up by simplifying. The built-in randomness should hopefully
    # take us somewhere new.
    tri.simplify()
    return


def _simplifyBasic(tri):
    """
    Uses 2-0 edge and 2-1 edge moves to monotonically reduce the number
    of tetrahedra in the given triangulation.

    This helper routine is used in randomise() to perform intermediate
    simplifications in a way that avoids undoing previously-performed random
    2-3 moves.

    If tri is currently oriented, then this routine guarantees to preserve
    the orientation.

    If no 2-0 or 2-1 moves are available, then the ambient triangulation
    will remain entirely untouched.

    Adapted from SnapPea's check_for_cancellation().

    Returns:
        True if and only if tri was successfully simplified.
    """
    changed = False     # Has anything changed ever?    (Return value.)
    changedNow = True   # Did we just change something? (Loop control.)
    while changedNow:
        changedNow = False
        for edge in tri.edges():
            # Try a 2-0 edge move.
            if tri.move20(edge):
                changedNow = True
                break

            # Try a 2-1 edge move.
            if tri.move21( edge, 0 ):
                changedNow = True
                break
            if tri.move21( edge, 1 ):
                changedNow = True
                break
        if changedNow:
            changed = True

    # Nothing further we can do.
    return changed
