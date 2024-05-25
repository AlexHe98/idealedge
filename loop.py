"""
Ideal loops for representing torus boundary components of a 3-manifold.
"""
from regina import *


def persistentLocation(face):
    """
    Returns a persistent identifier for the location of the given face.

    In detail, this routine returns a pair (t,f), where:
    --> t is the index of a tetrahedron meeting the given face; and
    --> f is a face number of this tetrahedron that corresponds to the given
        face.
    This identifier remains valid as long as tetrahedron indices remain
    unchanged (for example, as long as no tetrahedra are ever deleted).
    """
    emb = face.embedding(0)
    return ( emb.tetrahedron().index(), emb.face() )


class IdealLoop:
    """
    A sequence of edges representing an embedded ideal loop in a 3-manifold
    triangulation.
    """
    def __init__( self, edges ):
        """
        Creates an ideal loop from the given list of edges.

        Raises ValueError if the given list of edges does not form an
        embedded closed loop, or if the order of the edges in the given list
        does not match the order in which the edges appear in the loop.

        Pre-condition:
        --> The given list of edges is nonempty, and consists of edges that
            all belong to the same 3-manifold triangulation.
        """
        edge = edges[0]
        self._tri = edge.triangulation()
        self._edgeIndices = []
        self._vertIndices = set()
        self._tails = []
        firstVert = edge.vertex(0)
        lastVert = edge.vertex(0)
        error = False
        for edge in edges:
            self._edgeIndices.append( edge.index() )
            self._vertIndices.add( lastVert.index() )
            if edge.vertex(0) == lastVert:
                self._tails.append(0)
                lastVert = edge.vertex(1)
            elif edge.vertex(1) == lastVert:
                self._tails.append(1)
                lastVert = edge.vertex(0)
            else:
                error = True
                break
        if ( error or ( lastVert != firstVert ) or
                ( len( self._vertIndices ) != len(edges) ) ):
            msg = ( "Sequence of edges does not describe an embedded " +
                    "closed loop." )
            raise ValueError(msg)
        return

    def __len__(self):
        return len( self._edgeIndices )

    def __iter__(self):
        return iter( self._edgeIndices )

    def triangulation(self):
        """
        Returns the triangulation that contains this ideal loop.
        """
        return self._tri

    def weight( self, surf ):
        """
        Returns the number of times this ideal loop intersects the given
        normal surface surf.

        Pre-condition:
        --> The given normal surface is embedded in self.triangulation().
        """
        wt = 0
        for i in self._edgeIndices:
            wt += surf.edgeWeight(i).safeLongValue()
        return wt

    def components( self, surf ):
        """
        Returns a list describing the components into which the given normal
        surface surf splits this ideal loop.

        In detail, each item of the returned list is a list of edge segments.
        """
        lastComponent = []
        splitIndex = None
        for i in range( len(self) ):
            edgeIndex = self._edgeIndices[i]
            wt = surf.edgeWeight(edgeIndex).safeLongValue()
            if wt > 0:
                if self._tails[i] == 0:
                    lastComponent.append( ( edgeIndex, 0 ) )
                    headSeg = ( edgeIndex, wt )
                else:
                    lastComponent.append( ( edgeIndex, wt ) )
                    headSeg = ( edgeIndex, 0 )
                splitIndex = i
                break
            else:
                lastComponent.append( ( edgeIndex, 0 ) )
        if splitIndex is None:
            return [lastComponent]

        # The given surf splits this ideal loop into multiple components, so
        # we need to do a bit more work.
        components = []
        while splitIndex is not None:
            for seg in range( 1, wt ):
                components.append( [ ( edgeIndex, seg ) ] )
            nextComponent = [headSeg]
            continuation = splitIndex + 1
            splitIndex = None
            for i in range( continuation, len(self) ):
                edgeIndex = self._edgeIndices[i]
                wt = surf.edgeWeight(edgeIndex).safeLongValue()
                if wt > 0:
                    if self._tails[i] == 0:
                        nextComponent.append( ( edgeIndex, 0 ) )
                        headSeg = ( edgeIndex, wt )
                    else:
                        nextComponent.append( ( edgeIndex, wt ) )
                        headSeg = ( edgeIndex, 0 )
                    splitIndex = i
                    components.append(nextComponent)
                    break
                else:
                    nextComponent.append( ( edgeIndex, 0 ) )
        components.append( [ *nextComponent, *lastComponent ] )
        return components
