"""
Find the ideal edge after crushing a separating annulus in a one-vertex
triangulation of a 3-manifold with torus boundary.
"""
from regina import *


def _faceNumberings( emb, ind ):
    """
    Assuming that emb contains exactly two embeddings of a single edge e, and
    that ind is a 2-element list or tuple whose entries are both either 2 or
    3, returns the triangle indices and vertex numberings for each triangle
    that appears as "face ind[i] of emb[i]", for each i in {0,1}.

    To be precise, this routine returns a pair (f, v) such that:
    --> for i in {0,1}, f[i] is the index of the triangle given by
            emb[i].tetrahedron().triangle(
                emb[i].vertices()[ ind[i] ] );
    --> for i,j in {0,1}, v[i][j] is the vertex number of triangle f[i]
        corresponding to vertex j of the ede e; and
    --> for i in {0,1}, v[i][2] is the vertex number of triangle f[i] that
        is opposite the edge e.
    """
    f = []
    v = []
    for i in range(2):
        # Compute f[i].
        tet = emb[i].tetrahedron()
        ver = emb[i].vertices()
        face = tet.triangle( ver[ ind[i] ] )
        f.append( face.index() )

        # Compute v[i].
        #TODO Use tet.faceMapping()...
        if ind[i] == 2:
            p = Perm4(2,3)
        elif ind[i] == 3:
            p = Perm4()
        else:
            raise ValueError( "Entries of ind must be either 2 or 3." )
        faceEmb = TriangleEmbedding3( tet, ver*p )
        perm = faceEmb.vertices().inverse() * emb[i].vertices() * p[i]
        #TODO
        pass
    #TODO
    pass


def idealEdge(surf):
    """
    Assuming that surf is a separating annulus in a one-vertex triangulation
    of a 3-manifold with torus boundary, constructs the triangulation T given
    by crushing surf and returns the edge of T that corresponds to the ideal
    edge.
    
    Warning: This routine does *not* check that the input is sensible.
    """
    #TODO
    # Start by finding maximum-weight boundary edge.
    tri = surf.triangulation()
    bdry = tri.boundaryComponent(0)
    maxWt = 0
    edge = None
    for e in bdry.edges():
        wt = surf.edgeWeight( e.index() )
        if wt > maxWt:
            maxWt = wt
            edge = e

    #TODO
    #
    #TODO
    
    #TODO
    # Map tetrahedron indices to the new indices after crushing. 
    crushIndex = dict()
    #TODO
    
    #TODO
    #
    crushed = surf.crush()
    #TODO

    #TODO
    pass
