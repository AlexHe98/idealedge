"""
Speed up quad vertex normal surface enumeration.
"""
from enum import Enum, auto
from multiprocessing import Process, Pipe
from time import sleep
from regina import *
from loop import BoundsDisc
from triloops import EdgeIdealTriangulation, TriangulationWithBoundaryLoops
#TODO Consider unifying this implementation with the enumerations performed
#   as part of the knot factorisation algorithm.
#TODO The current implementation supports Triangulation3, but we never
#   actually use this because it requires dealing with dead ends, which makes
#   the interaction with the overall logic of the algorithm more subtle.


def findQuadVertexSurface( tri, identifyAcceptableSurface,
                          runParallelEnumerations=True ):
    """
    Searches for a quad vertex normal surface S in a triangulation T of the
    same 3-manifold as tri, such that identifyAcceptableSurface(T, S) is not
    None.

    If, for whatever reason, you require a surface specifically in the given
    triangulation (and not in some different triangulation of the same
    3-manifold), you should set runParallelEnumerations to False.

    The given tri must be an instance of either EdgeIdealTriangulation,
    TriangulationWithBoundaryLoops, or Regina's Triangulation3.

    If this routine finds an acceptable quad vertex surface S in a
    triangulation T of the same 3-manifold as tri, then it will return the
    triple (T, S, R), where R is the return value of
        identifyAcceptableSurface(T, S).
    Otherwise, this routine returns (T, None, None), where T is a
    triangulation of the same 3-manifold as tri in which it is certified that
    no acceptable quad vertex surface exists.

    If runParallelEnumerations is True (the default), then addition to
    running a main enumeration on tri (which is guaranteed to terminate, but
    which might be slow even if an acceptable surface exists), this routine
    will also start parallel processes which attempt to randomise tri in
    search of a different triangulation (of the same 3-manifold) for which an
    enumeration can terminate more quickly. In cases where an acceptable
    surfaces does indeed exist, this is often an effective method for
    speeding up the search.
    """
    if isinstance( tri, EdgeIdealTriangulation ):
        triType = _TriType.EDGE_IDEAL
        blueprint = tri.blueprint()
        underlyingTri = tri.triangulation()
    elif isinstance( tri, TriangulationWithBoundaryLoops ):
        triType = _TriType.BDRY_LOOPS
        blueprint = tri.blueprint()
        underlyingTri = tri.triangulation()
    elif isinstance( tri, Triangulation3 ):
        triType = _TriType.REGINA
        blueprint = tri.tightEncoding()
        underlyingTri = tri
    else:
        raise TypeError( "enumerationParallel() requires an instance of " +
                        "either EdgeIdealTriangulation or Triangulation3" )
    if runParallelEnumerations:
        # Set up a child process to repeatedly randomise the given triangulation,
        # and send the randomised triangulations to another child process that
        # runs alternate enumerations.
        randomiseReceiver, randomiseSender = Pipe(False)
        randomiseProcess = Process( target=_perpetualRandomise,
                args=( blueprint, randomiseSender ) )
        randomiseProcess.start()

        # Set up a child process to run the alternate enumerations.
        alternateReceiver, alternateSender = Pipe(False)
        alternateProcess = Process( target=_indefiniteEnumerate,
                args=( randomiseReceiver, alternateSender,
                      triType, identifyAcceptableSurface ) )
        alternateProcess.start()

    # Run the main enumeration.
    #NOTE As of Regina 7.4, NS_QUAD has been deprecated, and replaced with
    #   NormalCoords.Quad.
    enumeration = TreeEnumeration( underlyingTri, NormalCoords.Quad )
    while True:
        if runParallelEnumerations:
            # Has the randomiseProcess detected a loop which bounds a disc?
            if not randomiseProcess.is_alive():
                # Make sure to clean up child processes before raising
                # BoundsDisc.
                alternateProcess.terminate()
                randomiseProcess.join()
                alternateProcess.join()
                raise BoundsDisc()

            # Has the altenatveProcess given an answer?
            if alternateReceiver.poll():
                # Make sure to clean up child processes before returning the
                # answer from the alternateProcess.
                randomiseProcess.terminate()
                alternateProcess.join()
                randomiseProcess.join()
                #TODO For now, we ignore the number of attempts.
                newBlueprint, surfCoords, surfDesc, _ =\
                        alternateReceiver.recv()
                if triType == _TriType.EDGE_IDEAL:
                    newTri = EdgeIdealTriangulation.fromBlueprint(
                            *newBlueprint )
                    newUnderlyingTri = newTri.triangulation()
                elif triType == _TriType.BDRY_LOOPS:
                    newTri = TriangulationWithBoundaryLoops.fromBlueprint(
                            *newBlueprint )
                    newUnderlyingTri = newTri.triangulation()
                else:
                    newTri = Triangulation3.tightDecoding(newBlueprint)
                    newUnderlyingTri = newTri
                if surfCoords is None:
                    # No acceptable surface exists.
                    return ( newTri, None, None )
                foundSurf = NormalSurface(
                        newUnderlyingTri, NormalCoords.Standard, surfCoords )
                return ( newTri, foundSurf, surfDesc )

        # Continue with the main enumeration.
        if not enumeration.next():
            # No acceptable surfaces.
            ans = ( tri, None, None )
            break
        surf = enumeration.buildSurface()
        surfDesc = identifyAcceptableSurface( tri, surf )
        if surfDesc is not None:
            # Found an acceptable surface!
            ans = ( tri, surf, surfDesc )
            break

    # Clean up child processes before returning.
    if runParallelEnumerations:
        alternateProcess.terminate()
        randomiseProcess.terminate()
        alternateProcess.join()
        randomiseProcess.join()
    return ans


class _TriType(Enum):
    EDGE_IDEAL = auto()
    BDRY_LOOPS = auto()
    REGINA = auto()
    pass


def _perpetualRandomise( blueprint, sender ):
    attempts = 0
    if isinstance( blueprint, str ):
        tri = Triangulation3.tightDecoding(blueprint)
        size = tri.size()
        while True:
            attempts += 1
            randomise(tri)
            if tri.size() <= size:
                # Send randomised tri.
                sender.send( ( tri.tightEncoding(), attempts ) )
    else:
        edgeIdealTri = EdgeIdealTriangulation.fromBlueprint(*blueprint)
        size = edgeIdealTri.triangulation().size()
        while True:
            attempts += 1
            try:
                edgeIdealTri.randomise()    # Might raise BoundsDisc.
            except BoundsDisc:
                # This is the only route by which this routine can terminate
                # early, so callers can use early termination as a
                # certificate that some ideal loop bounds a disc.
                return
            if edgeIdealTri.triangulation().size() <= size:
                # Send randomised edgeIdealTri.
                sender.send( ( edgeIdealTri.blueprint(), attempts ) )
    return


def _indefiniteEnumerate(
        receiver, sender, triType, identifyAcceptableSurface ):
    searches = 0
    while not receiver.poll():
        sleep(0.01)
    blueprint, attempts = receiver.recv()
    if triType == _TriType.EDGE_IDEAL:
        tri = EdgeIdealTriangulation.fromBlueprint(*blueprint)
        underlyingTri = tri.triangulation()
    elif triType == _TriType.BDRY_LOOPS:
        tri = TriangulationWithBoundaryLoops.fromBlueprint(*blueprint)
        underlyingTri = tri.triangulation()
    else:
        tri = Triangulation3.tightDecoding(blueprint)
        underlyingTri = tri
    #NOTE As of Regina 7.4, NS_QUAD has been deprecated, and replaced with
    #   NormalCoords.Quad.
    enumeration = TreeEnumeration( underlyingTri, NormalCoords.Quad )
    while True:
        # From empirical observations, 20 searches seems to be a reasonable
        # cut-off for when it is probably worthwhile to give up and attempt
        # a new enumeration on a different triangulation instead.
        if searches > 20 and receiver.poll():
            searches = 0
            blueprint, attempts = receiver.recv()
            if triType == _TriType.EDGE_IDEAL:
                tri = EdgeIdealTriangulation.fromBlueprint(*blueprint)
                underlyingTri = tri.triangulation()
            elif triType == _TriType.BDRY_LOOPS:
                tri = TriangulationWithBoundaryLoops.fromBlueprint(*blueprint)
                underlyingTri = tri.triangulation()
            else:
                tri = Triangulation3.tightDecoding(blueprint)
                underlyingTri = tri
            #NOTE As of Regina 7.4, NS_QUAD has been deprecated, and replaced
            #   with NormalCoords.Quad.
            enumeration = TreeEnumeration( underlyingTri, NormalCoords.Quad )

        # Get the next surface and check whether it is acceptable.
        searches += 1
        if not enumeration.next():
            # No acceptable surface exists.
            sender.send( ( blueprint, None, None, attempts ) )
            return
        surf = enumeration.buildSurface()
        surfDesc = identifyAcceptableSurface( tri, surf )
        if surfDesc is not None:
            # Found an acceptable surface!
            sender.send( ( blueprint, _normalSurfacePythonVector(surf),
                          surfDesc, attempts ) )
            return


def _normalSurfacePythonVector(surf):
    return [ i.pythonValue() for i in surf.vector() ]
