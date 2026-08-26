"""
Recognition of bounded orientable Seifert fibred spaces.
"""


def recogniseSFS(tri):
    """
    Determines whether the given triangulation is a bounded orientable
    Seifert fibred space, and if so returns a Seifert fibration.

    If the triangulation is not a Seifert fibred space, then this routine
    returns None.

    This routine requires that tri is bounded and orientable, and raises
    ValueError if these conditions are not satisfied.

    Warning:
        The algorithms used in this routine rely on normal surface theory,
        and so might be very slow for larger triangulations (although faster
        tests are used where possible).
    """
    if not tri.isValid():
        # This rules out:
        #   --> Invalid edges
        #   --> Vertex links which are bounded and not discs
        raise ValueError(
                "recogniseSFS() requires a valid triangulation" )
    if not tri.isOrientable():
        raise ValueError(
                "recogniseSFS() requires an orientable triangulation" )
    if tri.isClosed() or tri.isIdeal():
        # This rules out:
        #   --> Triangulations in which all vertex links are 2-spheres
        #   --> Vertex links which are closed and not 2-spheres
        raise ValueError(
                "recogniseSFS() requires an orientable triangulation" )

    # At this point, we have an orientable 3-manifold with nonempty boundary,
    # where each boundary component is built entirely out of real boundary
    # triangles.
    for bc in tri.boundaryComponents():
        # Boundary components must all be tori to have a Seifert fibration.
        if ( bc.eulerChar() != 0 ) or ( not bc.isOrientable() ):
            return None
    tri.minimiseBoundary()

    # We now have a boundary-minimal triangulation of an orientable
    # 3-manifold whose boundary is a non-empty union of tori. This is where
    # the real work begins.

    #TODO Start with combinatorial recognition, and only fall back on normal
    #   surfaces and edge-ideal triangulations if this fails.
    #TODO Introduce a separate function which doesn't use combinatorial
    #   recognition, since this is what we want to experiment on.
    raise NotImplementedError()


class _SFSpaceRecognitionInvariants:
    """
    Internal invariants used by recogniseSFS() to help recover a complete
    description of a Seifert fibration.
    """
    def __init__(self):
        """
        Default initialises the invariants.

        Specifically, we will initially have the following:
        --> self.baseEuler() will be 0.
        --> self.fibres() will be an empty list.
        --> self.isBaseNonOrientable() will be False.
        """
        self._baseEuler = 0
        self._fibres = []
        self._isBaseNonOrbl = False
        return
    #TODO
