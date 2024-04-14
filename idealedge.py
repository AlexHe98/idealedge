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


def idealEdge(annulus):
    """
    Returns details of the ideal edge after crushing the given annulus.

    Specifically, this routine returns a pair (t,e), where t is the index
    (after crushing) of a tetrahedron that will meet the ideal edge after
    crushing, and e is an edge number of this tetrahedron that corresponds to
    the ideal edge.

    Warning:
    --> This routine assumes that the annulus is a separating surface in a
        one-vertex triangulation of a 3-manifold with torus boundary. This
        routine does not check that these conditions are satisfied.
    """
    tri = annulus.triangulation()

    # Split edges of tri along annulus.
    splitEdges = []
    for e in tri.edges():
        ei = e.index()
        wt = annulus.edgeWeight(ei)
        for s in range( wt + 1 ):
            splitEdges.append( (ei,s) )

    # Find "target" edges.
    targets = dict()
    for tet in tri.tetrahedra():
        teti = tet.index()
        hasQuads = False
        for q in range(3):
            if annulus.quads( teti, q ) > 0:
                hasQuads = True
                break
        if hasQuads:
            continue

        # No quads in tet, so there is a cell in the centre that survives
        # crushing. Find the edges of this cell that survive.
        for en in range(6):
            v = tet.edgeMapping(en)[0]
            ei = tet.edge(en).index()
            s = annulus.triangles( teti, v )
            targets[ (ei,s) ] = ( teti, en )

    # Starting edge for depth-first search.
    for e in tri.edges():
        ei = e.index()
        if ( not e.isBoundary() ) or ( annulus.edgeWeight(ei) <= 1 ):
            continue
        start = ( ei, 1 )

    # Find ideal edge using depth-first search.
    stack = [start]
    visited = set()
    while stack:
        ei, split = stack.pop()
        e = tri.edge(ei)
        for emb in e.embeddings():
            tet = emb.tetrahedron()
            teti = tet.index()
            en = emb.face()
            ver = emb.vertices()

            # Get normal coordinates incident to e in tet.
            f = [ annulus.triangles( teti, ver[i] ) for i in range(2) ]
            q = 0
            for qType in range(3):
                quads = annulus.quads( teti, qType )
                if quads > 0:
                    if qType not in { en, 5-en }:
                        q = quads
                    break

            # Does the current split edge belong to a parallel cell in this
            # tet?
            if split == f[0] or split == f[0] + q:
                continue

            # Is the parallel cell a triangular cell at vertex ver[0]?
            if split < f[0]:
                for enOther in range(6):
                    if enOther == en:
                        continue
                    verOther = tet.faceMapping(enOther)
                    if verOther[0] == ver[0]:
                        #TODO
                        pass
                    elif verOther[1] == ver[0]:
                        #TODO
                        pass
                #TODO
                continue
            #TODO

            # Is the parallel cell a quadrilateral cell?
            #TODO

            # The parallel cell must be a triangular cell at vertex ver[1].
            #TODO
            pass
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
