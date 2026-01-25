"""
Custom exceptions for EmbeddedLoop objects.
"""


class EmbeddedLoopException(Exception):
    pass


class NotLoop(EmbeddedLoopException):
    """
    Raised when attempting to build an embedded loop from a list of edges
    that does not describe a closed loop.
    """
    def __init__( self, edges ):
        indices = [ e.index() for e in edges ]
        msg = ( "The edge sequence {} does not describe ".format(indices) +
                "an embedded closed loop." )
        super().__init__(msg)
        return


class BoundsDisc(EmbeddedLoopException):
    """
    Raised when an embedded loop detects that it bounds an embedded disc.
    """
    def __init__(self):
        super().__init__( "Embedded loop bounds an embedded disc." )
        return
