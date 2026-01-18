"""
Perform 2-3, 3-2 and 2-0 moves while tracking how edges are relabelled.
"""
from regina import *
from edgelabel import EdgeLabelling


def edgeIndFromEmb(edgeEmbedding):
    """
    Returns the index of the underlying edge of the given EdgeEmbedding3
    object.
    """
    return edgeEmbedding.tetrahedron().edge( edgeEmbedding.edge() ).index()


def twoThree( triangle, edgeLab=None ):
    """
    Performs a 2-3 move about the given triangle, and returns an EdgeLabelling
    that tracks how edges were relabelled as a result of this move.

    More specifically, this routine tracks how edges are relabelled relative
    to the following "reference labelling" of (some or all of) the edges of
    triangle.triangulation():
    --> If the edgeLab parameter is omitted, then the default reference
        labelling assigns, to each edge e, the index e.index() and the
        embedding e.front().
    --> Otherwise, edgeLab must be an instance of EdgeLabelling with
        edgeLab.triangulation() == triangle.triangulation(), and the reference
        labelling is given by the (index, embedding) pairs that are specified
        by edgeLab.

    If the requested move is not legal, then triangle.triangulation() is left
    entirely untouched, and this routine returns None.

    Otherwise, if the move is legal, then this routine directly modifies the
    triangulation T := triangle.triangulation(). The returned EdgeLabelling r
    is structured as follows:
    --> For each index i in the reference labelling, which corresponds to some
        edge e in T before the move, r[i] will be an EdgeEmbedding3 object
        describing an embedding of e in T *after* performing the move. The
        embedding of e in r will have the same orientation as the embedding of
        e in the reference labelling.
    --> A 2-3 move also always creates a new edge, and r[i] will give an
        embedding of this new edge, where i is the largest negative index that
        is not already tracked by the reference labelling (typically, this
        will mean that i is -1).

    This routine will never modify edgeLab (if supplied).

    If triangle.triangulation() is currently oriented, then this orientation
    will be preserved by the requested 2-3 move.
    """
    tri = triangle.triangulation()
    if edgeLab is None:
        edgeLab = EdgeLabelling(tri)

    # Is the requested 2-3 move legal?
    #
    #NOTE Triangulation3.hasPachner(f) was introduced in Regina 7.4. In older
    #       versions of Regina, equivalent functionality (checking eligibility
    #       of the move, but not performing it) was provided by
    #       Triangulation3.pachner( f, True, False ).
    if not tri.hasPachner(triangle):
        return None

    # We need to work out the gluings that we need to perform before we make
    # any changes.
    doomed = []
    verts = []
    # Vertex numbering for tetrahedron doomed[i]. The vertex numbered
    # verts[i][3] forms an apex of the bipyramid.
    #
    #                        verts[i][3]
    #                             •
    #                            /|\
    #                           / | \
    #               verts[i][0]•--|--•verts[i][2]
    #                           \ | /
    #                            \|/
    #                             •
    #                        verts[i][1]
    #
    # The three new tetrahedra will be numbered survive+j, for j in {0,1,2}.
    # Vertices verts[0][j] and verts[0][3] of tetrahedron survive+j will each
    # be incident to an apex of the bipyramid. The faces opposite these two
    # vertices will form exposed faces of the bipyramid:
    #   --> face verts[0][j] of tetrahedron survive+j will correspond to face
    #       verts[0][j] of doomed[0]; and
    #   --> face verts[0][3] of tetrahedron survive+j will correspond to face
    #       verts[1][j] of doomed[1].
    #
    # Letting jj = (j+1)%3 and jm = (j-1)%3, tetrahedra survive+j and
    # survive+jj appear adjacent to each other in the bipyramid as follows:
    #
    #       Tetrahedron survive+jj           Tetrahedron survive+j
    #                                                         
    #            verts[0][3]                      verts[0][3]
    #                 •                                •
    #                /|\                              /|\
    #               / | \                            / | \
    #   verts[0][j]•--|--•verts[0][jm]  verts[0][jm]•--|--•verts[0][jj]
    #               \ | /                            \ | /
    #                \|/                              \|/
    #                 •                                •
    #            verts[0][jj]                     verts[0][j]
    #
    # That is, face verts[0][jj] of tetrahedron survive+j should be glued to
    # tetrahedron survive+jj by the permutation that swaps verts[0][j] and
    # verts[0][jj].
    for emb in triangle.embeddings():
        doomed.append( emb.simplex() )
        verts.append( emb.vertices() )
    newIndex = []   # New tetrahedron indices after performing the move.
    doomedIndices = sorted( [ d.index() for d in doomed ] )
    for k in range( tri.size() ):
        if k < doomedIndices[0]:
            newIndex.append(k)
        elif k < doomedIndices[1]:
            newIndex.append(k-1)
        else:
            newIndex.append(k-2)
    survive = tri.size() - 2
    gluings = []
    triGlu = doomed[0].adjacentGluing( verts[0][3] )
    for j in range(3):
        jj = (j+1)%3

        # Glue the new tetrahedron that will eventually be numbered survive+j
        # to the new tetrahedron that will eventually be numbered survive+jj.
        gluings.append(
                ( survive+j, verts[0][jj], survive+jj,
                    Perm4( verts[0][j], verts[0][jj] ) ) )

        # Glue triangle verts[0][j] of survive+j.
        oldAdj = doomed[0].adjacentSimplex( verts[0][j] )
        if oldAdj is not None:
            oldGlu = doomed[0].adjacentGluing( verts[0][j] )
            try:
                i = doomed.index(oldAdj)
            except ValueError:
                # We are gluing to a tetrahedron that survives.
                newAdj = newIndex[ oldAdj.index() ]
                gluings.append(
                        ( survive+j, verts[0][j], newAdj, oldGlu ) )
            else:
                # We are gluing to a doomed tetrahedron, which means that we
                # are actually supposed to glue to one of the new tetrahedra.
                jjj = verts[i].inverse()[ oldGlu[ verts[0][j] ] ]

                # Do the gluing if we haven't already done it from the other
                # side.
                if jjj >= j:
                    if i == 0:
                        gluings.append(
                                ( survive+j, verts[0][j],
                                    survive+jjj, oldGlu ) )
                    else:
                        newGlu = ( triGlu.inverse() *
                                Perm4( verts[1][3], verts[1][jjj] ) *
                                oldGlu )
                        gluings.append(
                                ( survive+j, verts[0][j],
                                    survive+jjj, newGlu ) )

        # Glue triangle verts[0][3] of survive+j.
        oldAdj = doomed[1].adjacentSimplex( verts[1][j] )
        if oldAdj is not None:
            oldGlu = doomed[1].adjacentGluing( verts[1][j] )
            try:
                i = doomed.index(oldAdj)
            except ValueError:
                # We are gluing to a tetrahedron that survives.
                newAdj = newIndex[ oldAdj.index() ]
                newGlu = oldGlu * Perm4( verts[1][3], verts[1][j] ) * triGlu
                gluings.append(
                        ( survive+j, verts[0][3], newAdj, newGlu ) )
            else:
                # We are gluing to a doomed tetrahedron, which means that we
                # are actually supposed to glue to one of the new tetrahedra.
                jjj = verts[i].inverse()[ oldGlu[ verts[1][j] ] ]

                # Do the gluing if we haven't already done it from the other
                # side.
                if jjj > j:
                    if i == 0:
                        newGlu = ( oldGlu *
                                Perm4( verts[1][3], verts[1][j] ) *
                                triGlu )
                        gluings.append(
                                ( survive+j, verts[0][3],
                                    survive+jjj, newGlu ) )
                    else:
                        newGlu = ( triGlu.inverse() *
                                Perm4( verts[1][3], verts[1][jjj] ) *
                                oldGlu *
                                Perm4( verts[1][3], verts[1][j] ) *
                                triGlu )
                        gluings.append(
                                ( survive+j, verts[0][3],
                                    survive+jjj, newGlu ) )

    # For each edge e that is tracked by the reference labelling, find one
    # tetrahedron that will meet e after we have performed the requested 2-3
    # move.
    newEdgeLocations = dict()
    for edgeInd in edgeLab:
        emb = edgeLab[edgeInd]
        oldTet = emb.tetrahedron()

        # For subsequent comments, let
        #   e := oldTet.edge( emb.edge() ).
        try:
            i = doomed.index(oldTet)
        except ValueError:
            # The oldTet survives.
            newEdgeLocations[edgeInd] = ( newIndex[ oldTet.index() ],
                                         emb.vertices() )
        else:
            # The tetrahedron oldTet is doomed, but this means that e will
            # meet one of the new tetrahedra.
            #
            # We now have two permutations of the vertices of oldTet:
            #   --> verts[i], which comes from the embedding of the input
            #       triangle in oldTet
            #   --> emb.vertices(), which comes from the embedding of the edge
            #       e in oldTet
            #
            # Also, recall that the new tetrahedra will be indexed by
            # survive+ii, ii in {0,1,2}, with vertices labelled as follows (in
            # the diagram, the front vertical edge is the new internal edge
            # introduced by the requested 2-3 move):
            #
            #                    Tetrahedron survive+ii
            #                                     
            #                         verts[0][3]
            #                              •
            #                             /|\
            #                            / | \
            #         verts[0][(ii-1)%3]•--|--•verts[0][(ii+1)%3]
            #                            \ | /
            #                             \|/
            #                              •
            #                         verts[0][ii]
            #
            p = verts[i].inverse() * emb.vertices()
            # By construction, verts[i][ p[j] ] == emb.vertices()[j].

            # Find suitable ii such that e will be incident to tetrahedron
            # survive+ii after performing the requested 2-3 move. This will
            # depend on whether e has one of its endpoints at an apex of the
            # bipyramid.
            apex = p.inverse()[3]
            if apex in {0,1}:
                # One of the endpoints of e is incident to an apex of the
                # bipyramid.
                j = p[ 1 - apex ]

                # We have e == oldTet.edge( verts[i][j], verts[i][3] ). In
                # particular, this means that after the 2-3 move, e will be
                # incident to the two new tetrahedra at the following indices:
                #   --> survive + ( (j+1) % 3 )
                #   --> survive + ( (j-1) % 3 )
                ii = (j+1) % 3
            else:   # apex in {2,3}
                # We have e == oldTet.edge( verts[i][p[0]], verts[i][p[1]] ).
                ii = p[ 5 - apex ]

            # Convert emb.vertices() into a vertex permutation that
            # describes how e will be embedded in tetrahedron survive+ii
            # after the requested 2-3 move.
            #
            # With ii chosen as above, this is the same regardless of the
            # value of apex.
            if i == 0:
                newVertPerm = emb.vertices()
            else:   # i == 1
                newVertPerm = ( Perm4( verts[0][3], verts[0][ii] ) *
                               triGlu.inverse() * emb.vertices() )

            newEdgeLocations[edgeInd] = ( survive+ii, newVertPerm )

    # Remove the two doomed tetrahedra, add three new tetrahedra, and then
    # perform the gluings that we just computed.
    for d in doomed:
        tri.removeTetrahedron(d)
    tri.newTetrahedra(3)
    for me, face, you, glu in gluings:
        tri.tetrahedron(me).join( face, tri.tetrahedron(you), glu )

    # How did the edges get relabelled?
    newLab = dict()
    for edgeInd in newEdgeLocations:
        newTetInd, newVertPerm = newEdgeLocations[edgeInd]
        newLab[edgeInd] = EdgeEmbedding3(
                tri.tetrahedron(newTetInd), newVertPerm )

    # Don't forget that we created a new edge! For each ii in {0,1,2}, the new
    # edge is embedded in tetrahedron survive+ii so that its endpoints lie on
    # the vertices numbered verts[0][3] and verts[0][ii].
    newRefTet = tri.tetrahedron(survive)
    newEdgeNum = Edge3.faceNumber( verts[0] * Perm4(1,3) )
    newEdgeInd = -1 + min( 0, *edgeLab.trackedIndices() )
    newLab[newEdgeInd] = EdgeEmbedding3(
            newRefTet, newRefTet.edgeMapping(newEdgeNum) )
    return EdgeLabelling( tri, newLab )


def threeTwo( edge, edgeLab=None ):
    """
    Performs a 3-2 move about the given edge, and returns an EdgeLabelling
    that tracks how edges were relabelled as a result of this move.

    More specifically, this routine tracks how edges are relabelled relative
    to the following "reference labelling" of (some or all of) the edges of
    edge.triangulation():
    --> If the edgeLab parameter is omitted, then the default reference
        labelling assigns, to each edge e, the index e.index() and the
        embedding e.front().
    --> Otherwise, edgeLab must be an instance of EdgeLabelling with
        edgeLab.triangulation() == edge.triangulation(), and the reference
        labelling is given by the (index, embedding) pairs that are specified
        by edgeLab.

    If the requested move is not legal, then edge.triangulation() is left
    entirely untouched, and this routine returns None.

    Otherwise, if the move is legal, then this routine directly modifies the
    triangulation T := edge.triangulation(). The returned EdgeLabelling r is
    structured as follows:
    --> The requested 3-2 move always destroys the given edge, so for any
        index i in the reference labelling that corresponds to the given edge,
        r[i] will be None.
    --> For every other index i in the reference labelling, which corresponds
        to some edge e in T (other than the given edge) before the move, r[i]
        will be an EdgeEmbedding3 object describing an embedding of e in T
        *after* performing the move. The embedding of e in r will have the
        same orientation as the embedding of e in the reference labelling.

    This routine will never modify edgeLab (if supplied).

    If edge.triangulation() is currently oriented, then this orientation will
    be preserved by the requested 3-2 move.
    """
    tri = edge.triangulation()
    if edgeLab is None:
        edgeLab = EdgeLabelling(tri)

    # Is the requested 3-2 move legal?
    #
    #NOTE Triangulation3.hasPachner(e) was introduced in Regina 7.4. In older
    #       versions of Regina, equivalent functionality (checking eligibility
    #       of the move, but not performing it) was provided by
    #       Triangulation3.pachner( e, True, False ).
    if not tri.hasPachner(edge):
        return None

    # We need to work out the gluings that we need to perform before we make
    # any changes.
    doomed = []
    verts = []
    for emb in edge.embeddings():
        doomed.append( emb.simplex() )
        verts.append( emb.vertices() )
    newIndex = [] # Value of newIndex[ doomed[i].index() ] is meaningless.
    doomedIndices = sorted( [ d.index() for d in doomed ] )
    for k in range( tri.size() ):
        if k < doomedIndices[0]:
            newIndex.append(k)
        elif k < doomedIndices[1]:
            newIndex.append(k-1)
        elif k < doomedIndices[2]:
            newIndex.append(k-2)
        else:
            newIndex.append(k-3)
    survive = tri.size() - 3
    gluings = [ ( survive, verts[1][1], survive+1,
        Perm4( verts[1][0], verts[1][1] ) ) ]
    # The two new tetrahedra will end up being numbered survive+j, for j in
    # {0,1}. These tetrahedra will be labelled in such a way that triangle
    # verts[1][j] of survive+j corresponds precisely to triangle verts[1][j]
    # of doomed[1], with precisely the same vertex numberings. Since this
    # doesn't tell us how doomed[1] meets the other two doomed tetrahedra,
    # we record this information now in the variable triGlu.
    #
    # With the vertices of doomed[1] labelled as in the diagram below, so that
    # the input edge is the front vertical edge, we will have the following:
    #   --> The new tetrahedron survive+0 will have triangle verts[1][0] in
    #       common with doomed[1], and will therefore appear as the lower
    #       tetrahedron in the bipyramid.
    #   --> The new tetrahedron survive+1 will have triangle verts[1][1] in
    #       common with doomed[1], and will therefore appear as the upper
    #       tetrahedron in the bipyramid.
    #
    #                            verts[1][0]
    #                                 •
    #                                /|\
    #                               / | \
    #                   verts[1][3]•--|--•verts[1][2]
    #                               \ | /
    #                                \|/
    #                                 •
    #                            verts[1][1]
    #
    triGlu = [
            doomed[1].adjacentGluing( verts[1][3] ),  # Gluing to doomed[0]
            None,
            doomed[1].adjacentGluing( verts[1][2] ) ] # Gluing to doomed[2]
    for j in range(2):
        # We have already glued triangle verts[1][1-j] of survive+j (to
        # triangle verts[1][j] of survive+(1-j)). We now glue the other
        # triangles of survive+j.

        # Glue triangle verts[1][j] of survive+j.
        oldAdj = doomed[1].adjacentSimplex( verts[1][j] )
        if oldAdj is not None:
            oldGlu = doomed[1].adjacentGluing( verts[1][j] )
            try:
                i = doomed.index(oldAdj)
            except ValueError:
                # We are gluing to a tetrahedron that survives.
                newAdj = newIndex[ oldAdj.index() ]
                gluings.append(
                        ( survive+j, verts[1][j], newAdj, oldGlu ) )
            else:
                # We are gluing to a doomed tetrahedron, which means that we
                # are actually supposed to glue to one of the new tetrahedra.
                jj = verts[i].inverse()[ oldGlu[ verts[1][j] ] ]

                # Do the gluing if we haven't already done it from the other
                # side.
                if jj >= j:
                    if i == 1:
                        # In this case, we know that jj == 1 and j == 0.
                        gluings.append(
                                ( survive, verts[1][0], survive+1, oldGlu ) )
                    else:
                        newGlu = ( triGlu[i].inverse() *
                                Perm4( verts[i][2+i//2], verts[i][jj] ) *
                                oldGlu )
                        gluings.append(
                                ( survive+j, verts[1][j],
                                    survive+jj, newGlu ) )

        # Glue triangle verts[1][2] of survive+j.
        oldAdj = doomed[2].adjacentSimplex( verts[2][j] )
        if oldAdj is not None:
            oldGlu = doomed[2].adjacentGluing( verts[2][j] )
            try:
                i = doomed.index(oldAdj)
            except ValueError:
                # We are gluing to a tetrahedron that survives.
                newAdj = newIndex[ oldAdj.index() ]
                newGlu = ( oldGlu *
                        Perm4( verts[2][3], verts[2][j] ) *
                        triGlu[2] )
                gluings.append(
                        ( survive+j, verts[1][2], newAdj, newGlu ) )
            else:
                # We are gluing to a doomed tetrahedron, which means that we
                # are actually supposed to glue to one of the new tetrahedra.
                jj = verts[i].inverse()[ oldGlu[ verts[2][j] ] ]

                # Do the gluing if we haven't already done it from the other
                # side.
                if jj >= j:
                    if i == 1:
                        # In this case, if jj == j then we would already have
                        # glued from the other side.
                        if jj > j:
                            # In this case, we know that jj == 1 and j == 0.
                            newGlu = ( oldGlu *
                                    Perm4( verts[2][3], verts[2][0] ) *
                                    triGlu[2] )
                            gluings.append(
                                    ( survive, verts[1][2],
                                        survive+1, newGlu ) )
                    else:
                        newGlu = ( triGlu[i].inverse() *
                                Perm4( verts[i][2+i//2], verts[i][jj] ) *
                                oldGlu *
                                Perm4( verts[2][3], verts[2][j] ) *
                                triGlu[2] )
                        gluings.append(
                                ( survive+j, verts[1][2],
                                    survive+jj, newGlu ) )

        # Glue triangle verts[1][3] of survive+j.
        oldAdj = doomed[0].adjacentSimplex( verts[0][j] )
        if oldAdj is not None:
            oldGlu = doomed[0].adjacentGluing( verts[0][j] )
            try:
                i = doomed.index(oldAdj)
            except ValueError:
                # We are gluing to a tetrahedron that survives.
                newAdj = newIndex[ oldAdj.index() ]
                newGlu = ( oldGlu *
                        Perm4( verts[0][2], verts[0][j] ) *
                        triGlu[0] )
                gluings.append(
                        ( survive+j, verts[1][3], newAdj, newGlu ) )
            else:
                # We are gluing to a doomed tetrahedron, which means that we
                # are actually supposed to glue to one of the new tetrahedra.
                jj = verts[i].inverse()[ oldGlu[ verts[0][j] ] ]

                # Do the gluing if we haven't already done it from the other
                # side.
                if jj > j:
                    # We know that jj == 1 and j == 0.
                    if i == 1:
                        newGlu = ( oldGlu *
                                Perm4( verts[0][2], verts[0][0] ) *
                                triGlu[0] )
                        gluings.append(
                                ( survive, verts[1][3], survive+1, newGlu ) )
                    else:
                        newGlu = ( triGlu[i].inverse() *
                                Perm4( verts[i][2+i//2], verts[i][1] ) *
                                oldGlu *
                                Perm4( verts[0][2], verts[0][0] ) *
                                triGlu[0] )
                        gluings.append(
                                ( survive, verts[1][3], survive+1, newGlu ) )

    # For each edge e (other than the input edge) that is tracked by the
    # reference labelling, find one tetrahedron that will meet e after we have
    # performed the requested 3-2 move.
    newEdgeLocations = dict()
    for edgeInd in edgeLab:
        emb = edgeLab[edgeInd]
        oldTet = emb.tetrahedron()
        if oldTet.edge( emb.edge() ) == edge:
            # This edge is destroyed by the 3-2 move. Ignoring it now will
            # mean that the returned EdgeLabelling will automatically stop
            # tracking this edge.
            continue

        # For subsequent comments, let
        #   e := oldTet.edge( emb.edge() ).
        try:
            i = doomed.index(oldTet)
        except ValueError:
            # The oldTet survives.
            newEdgeLocations[edgeInd] = ( newIndex[ oldTet.index() ],
                                         emb.vertices() )
        else:
            # The tetrahedron oldTet is doomed, but this means that e will
            # meet one of the new tetrahedra.
            #
            # We now have two permutations of the vertices of oldTet:
            #   --> verts[i], which comes from the embedding of the input edge
            #       in oldTet
            #   --> emb.vertices(), which comes from the embedding of the edge
            #       e in oldTet
            p = verts[i].inverse() * emb.vertices()
            # By construction, verts[i][ p[j] ] == emb.vertices()[j].

            # Find suitable ii such that e will be incident to tetrahedron
            # survive+ii after performing the requested 3-2 move. This will
            # depend on whether e has one of its endpoints at an apex of the
            # bipyramid.
            top = p.inverse()[0]
            bot = p.inverse()[1]
            # Note that at most one of top or bot is in {0,1} (otherwise we
            # would have e == edge).
            if top in {0,1}:
                # The edge e is incident to the apex at vertex 0 of the input
                # edge.
                ii = 1
            elif bot in {0,1}:
                # The edge e is incident to the apex at vertex 1 of the input
                # edge.
                ii = 0
            else:
                # The edge e runs along the "equator" of the bipyramid, so
                # either choice of ii in {0,1} would work.
                ii = 0

            # Convert emb.vertices() into a vertex permutation that
            # describes how e will be embedded in tetrahedron survive+ii
            # after the requested 2-3 move.
            #
            # With ii chosen as above, this is the same regardless of the
            # values of top and bot.
            #
            # Recall that triangle verts[1][ii] of tetrahedron survive+ii
            # corresponds precisely to triangle verts[1][ii] of doomed[1],
            # with precisely the same vertex numberings.
            if i == 1:
                newVertPerm = emb.vertices()
            elif i == 0:
                #
                #     Tetrahedron doomed[1]          Tetrahedron doomed[0]
                #
                #          verts[1][ii]                   verts[0][ii]
                #               •                              •
                #              /|\                            /|\
                #             / | \                          / | \
                # verts[1][3]•--|--•verts[1][2]  verts[0][3]•--|--•verts[0][2]
                #             \ | /                          \ | /
                #              \|/                            \|/
                #               •                              •
                #          verts[1][1-ii]                 verts[0][1-ii]
                #
                newVertPerm = ( Perm4( verts[1][ii], verts[1][3] ) *
                               triGlu[0].inverse() * emb.vertices() )
            else:   # i == 2
                #
                #     Tetrahedron doomed[2]          Tetrahedron doomed[1]
                #
                #          verts[2][ii]                   verts[1][ii]
                #               •                              •
                #              /|\                            /|\
                #             / | \                          / | \
                # verts[2][3]•--|--•verts[2][2]  verts[1][3]•--|--•verts[1][2]
                #             \ | /                          \ | /
                #              \|/                            \|/
                #               •                              •
                #          verts[2][1-ii]                 verts[1][1-ii]
                #
                newVertPerm = ( Perm4( verts[1][ii], verts[1][2] ) *
                               triGlu[2].inverse() * emb.vertices() )

            newEdgeLocations[edgeInd] = ( survive+ii, newVertPerm )

    # Remove the three doomed tetrahedra, add two new tetrahedra, and then
    # perform the gluings that we just computed.
    removeIndex = edge.index() # Need to remember this for later.
    for d in doomed:
        tri.removeTetrahedron(d)
    tri.newTetrahedra(2)
    for me, face, you, glu in gluings:
        tri.tetrahedron(me).join( face, tri.tetrahedron(you), glu )

    # How did the edges get relabelled?
    newLab = dict()
    for edgeInd in newEdgeLocations:
        newTetInd, newVertPerm = newEdgeLocations[edgeInd]
        newLab[edgeInd] = EdgeEmbedding3(
                tri.tetrahedron(newTetInd), newVertPerm )

    # All done!
    return EdgeLabelling( tri, newLab )


def twoZero( edge, edgeLab=None ):
    """
    Performs a 2-0 move about the given edge, and returns an EdgeLabelling
    that tracks how edges were relabelled as a result of this move.

    More specifically, this routine tracks how edges are relabelled relative
    to the following "reference labelling" of (some or all of) the edges of
    edge.triangulation():
    --> If the edgeLab parameter is omitted, then the default reference
        labelling assigns, to each edge e, the index e.index() and the
        embedding e.front().
    --> Otherwise, edgeLab must be an instance of EdgeLabelling with
        edgeLab.triangulation() == edge.triangulation(), and the reference
        labelling is given by the (index, embedding) pairs that are specified
        by edgeLab.

    If the requested move is not legal, then edge.triangulation() is left
    entirely untouched, and this routine returns None.

    Otherwise, if the move is legal, then this routine directly modifies the
    triangulation T := edge.triangulation(). The returned EdgeLabelling r is
    structured as follows:
    --> The requested 2-0 move always destroys the given edge, so for any
        index i in the reference labelling that corresponds to the given edge,
        r[i] will be None.
    --> For every other index i in the reference labelling, which corresponds
        to some edge e in T (other than the given edge) before the move, r[i]
        will be an EdgeEmbedding3 object describing an embedding of e in T
        *after* performing the move. The embedding of e in r will have the
        same orientation as the embedding of e in the reference labelling.

    This routine will never modify edgeLab (if supplied).

    Note also that a 2-0 move merges two edges into a single new edge e. If
    the reference labelling tracks indices, say i and j, for both of the
    merged edges, then in the returned EdgeLabelling r, we will have that r[i]
    and r[j] give two (possibly equal, possibly distinct) EdgeEmbedding3
    objects corresponding to the same new edge e.

    If edge.triangulation() is currently oriented, then this orientation will
    be preserved by the requested 2-0 move.
    """
    tri = edge.triangulation()
    if edgeLab is None:
        edgeLab = EdgeLabelling(tri)

    # Is the requested 2-0 move legal?
    #
    #NOTE Triangulation3.has20(e) was introduced in Regina 7.4. In older
    #       versions of Regina, equivalent functionality (checking eligibility
    #       of the move, but not performing it) was provided by
    #       Triangulation3.twoZeroMove( e, True, False ).
    if not tri.has20(edge):
        return None

    # How will the tetrahedra in tri get renumbered after we perform the
    # requested 2-0 move?
    doomedIndices = sorted(
            [ emb.simplex().index() for emb in edge.embeddings() ] )
    newIndex = [] # Value of newIndex[ doomedIndices[i] ] is meaningless.
    for k in range( tri.size() ):
        if k < doomedIndices[0]:
            newIndex.append(k)
        elif k < doomedIndices[1]:
            newIndex.append(k-1)
        else:
            newIndex.append(k-2)

    # For each edge e (other than the input edge) that is tracked by the
    # reference labelling, find one tetrahedron that will meet e after we have
    # performed the requested 2-0 move.
    newEdgeLocations = dict()
    for edgeInd in edgeLab:
        emb = edgeLab[edgeInd]
        tet = emb.tetrahedron()
        if tet.edge( emb.edge() ) == edge:
            # This edge is destroyed by the 2-0 move. Ignoring it now will
            # mean that the returned EdgeLabelling will automatically stop
            # tracking this edge.
            continue

        # Go through the embeddings of e := tet.edge( emb.edge() ), and
        # look for a tetrahedron that survives the 2-0 move.
        found = False
        for otherEmb in tet.edge( emb.edge() ).embeddings():
            oldTetInd  = otherEmb.tetrahedron().index()
            if oldTetInd in doomedIndices:
                continue

            # This one survives!
            #
            # If otherEmb.vertices() induces the opposite orientation on e to
            # the one given by emb.vertices(), then we will need to flip the
            # orientation.
            found = True
            otherVer = otherEmb.vertices()
            if otherVer[0] == emb.vertices()[0]:
                newVertPerm = otherVer
            else:
                newVertPerm = otherVer * Perm4(1,0,3,2)
            newEdgeLocations[edgeInd] = ( newIndex[oldTetInd], newVertPerm )
            break
        if found:
            continue

        # There is only one way we can fail to find a surviving tetrahedron.
        # To see this, recall that the two tetrahedra involved in the 2-0 move
        # form a square pillow, and we have already checked that the move is
        # legal. Since e is not incident to a surviving tetrahedron, it must
        # be the case that e is one of the two edges "opposite" the input
        # edge; moreover, we must have folded the two adjacent faces across e
        # so that e forms an internal edge of degree one. The 2-0 move will
        # merge e with the other opposite edge, and so we can simply go
        # through the tetrahedra incident to this other edge and thereby find
        # a suitable surviving tetrahedron.
        #
        # We first find the other opposite edge by examining the two
        # tetrahedra that are involved in the 2-0 move.
        #
        #           tet                         otherOppTet
        #
        #          ver[2]                      otherOppVer[2]
        #            •                               •
        #           /|\                             /|\
        #          / | \                           / | \
        #   ver[0]•--|--•ver[1]     otherOppVer[1]•--|--•otherOppVer[0]
        #          \ | /                           \ | /
        #           \|/                             \|/
        #            •                               •
        #          ver[3]                      otherOppVer[3]
        #
        ver = emb.vertices()
        gluing = tet.adjacentGluing( ver[0] )
        otherOppTet = tet.adjacentTetrahedron( ver[1] )
        otherOppVer = gluing * ver
        otherOppEdgeNum = Edge3.faceNumber(otherOppVer)
        otherOppEdge = otherOppTet.edge(otherOppEdgeNum)

        # Find a surviving tetrahedron incident to the otherOppEdge.
        fixOrientation = None
        surviveEmb = None
        for candidateEmb in otherOppEdge.embeddings():
            candidateTet = candidateEmb.tetrahedron()
            if candidateTet == otherOppTet:
                # Implies candidateTet.index() in doomedIndices.
                if candidateEmb.edge() != otherOppEdgeNum:
                    continue

                # Compare the vertex permutation given by candidateEmb with
                # the permutation otherOppVer. This will tell us whether the
                # 2-0 move will merge e and otherOppEdge with the same or
                # opposite orientations.
                if candidateEmb.vertices()[0] == otherOppVer[0]:
                    fixOrientation = Perm4()
                else:
                    fixOrientation = Perm4(1,0,3,2)

                # If we have already determined surviveEmb, then we are ready
                # to determine the newEdgeLocations data for e.
                if surviveEmb is not None:
                    break
            elif candidateTet.index() not in doomedIndices:
                # This one survives!
                if surviveEmb is None:
                    surviveEmb = candidateEmb

                # If we have already determined fixOrientation, then we are
                # ready to determine the newEdgeLocations data for e.
                if fixOrientation is not None:
                    break

        newEdgeLocations[edgeInd] = (
                newIndex[ surviveEmb.tetrahedron().index() ],
                surviveEmb.vertices() * fixOrientation )

    # Perform the 2-0 move, and work out how the edges were relabelled.
    #
    #NOTE Triangulation3.move20(e) was introduced in Regina 7.4. In older
    #       versions of Regina, 2-0 moves were performed using
    #           Triangulation3.twoZeroMove(e).
    tri.move20(edge)
    newLab = dict()
    for edgeInd in newEdgeLocations:
        newTetInd, newVertPerm = newEdgeLocations[edgeInd]
        newLab[edgeInd] = EdgeEmbedding3(
                tri.tetrahedron(newTetInd), newVertPerm )

    # All done!
    return EdgeLabelling( tri, newLab )


def fourFour( edge, newAxis, edgeLab=None ):
    """
    Performs a 4-4 move about the given edge, and returns an EdgeLabelling
    that tracks how edges were relabelled as a result of this move.

    More specifically, this routine tracks how edges are relabelled relative
    to the following "reference labelling" of (some or all of) the edges of
    edge.triangulation():
    --> If the edgeLab parameter is omitted, then the default reference
        labelling assigns, to each edge e, the index e.index() and the
        embedding e.front().
    --> Otherwise, edgeLab must be an instance of EdgeLabelling with
        edgeLab.triangulation() == edge.triangulation(), and the reference
        labelling is given by the (index, embedding) pairs that are specified
        by edgeLab.

    If the 4-4 move is legal, then the given edge forms the axis of a
    four-tetrahedron octahedron. There are two other choices of axis for this
    octahedron, and the requested 4-4 move retriangulates this octahedron to
    use one of these other choices of axis edge. Specifically, number the
    original four tetrahedra 0, 1, 2 and 3 according to the order described
    by e.embeddings(). If newAxis is 0, then the new axis will separate
    tetrahedra 0 and 1 from tetrahedra 2 and 3. If newAxis is 1, then the new
    axis will separate tetrahedra 1 and 2 from tetrahedra 3 and 0.

    If the requested move is not legal, then edge.triangulation() is left
    entirely untouched, and this routine returns None.

    Otherwise, if the move is legal, then this routine directly modifies the
    triangulation T := edge.triangulation(). The returned EdgeLabelling r is
    structured as follows:
    --> The requested 4-4 move always destroys the given edge, so for any
        index i in the reference labelling that corresponds to the given edge,
        r[i] will be None.
    --> For every other index i in the reference labelling, which corresponds
        to some edge e in T (other than the given edge) before the move, r[i]
        will be an EdgeEmbedding3 object describing an embedding of e in T
        *after* performing the move. The embedding of e in r will have the
        same orientation as the embedding of e in the reference labelling.
    --> As mentioned above, a 4-4 move also always creates a new axis edge,
        and r[i] will give an embedding of this new edge, where i is the
        largest negative index that is not already tracked by the reference
        labelling (typically, this will mean that i is -1).

    This routine will never modify edgeLab (if supplied).

    If edge.triangulation() is currently oriented, then this orientation will
    be preserved by the requested 4-4 move.
    """
    tri = edge.triangulation()
    if edgeLab is None:
        edgeLab = EdgeLabelling(tri)

    # Is the requested 4-4 move legal?
    #
    #NOTE Triangulation3.has44( e, ax ) was introduced in Regina 7.4. In older
    #       versions of Regina, equivalent functionality (checking eligibility
    #       of the move, but not performing it) was provided by
    #       Triangulation3.fourFourMove( e, ax, True, False ).
    if not tri.has44( edge, newAxis ):
        return None

    # Find the doomed tetrahedra.
    doomed = []
    verts = []
    for emb in edge.embeddings():
        doomed.append( emb.simplex() )
        verts.append( emb.vertices() )

    # The requested 4-4 move is a composition of the following two moves.
    #   --> A 2-3 move on a triangle f23 determined by the choice of newAxis:
    #       --- If newAxis == 0, then we choose f23 to be the triangle between
    #           doomed[0] and doomed[1].
    #       --- If newAxis == 1, then we choose f23 to be the triangle between
    #           doomed[1] and doomed[2].
    #   --> A 3-2 move on the input edge.
    # Note that doomed[3] is never involved in the initial 2-3 move, so we can
    # keep track of the input edge via its embedding in doomed[3].
    if newAxis == 0:
        f23 = doomed[0].triangle( verts[0][2] )
    else:
        f23 = doomed[1].triangle( verts[1][2] )
    e32 = edge.embedding(3).edge()
    relab = twoThree(f23)
    return threeTwo( doomed[3].edge(e32), relab )


#TODO Reimplement using EdgeLabelling.
def twoOne( edge, edgeEnd ):
    """
    Performs a 2-1 move about the given edge, and returns a dictionary r that
    describes how the edges were renumbered.

    If the 2-1 move is legal, then the given edge forms the internal edge of
    a snapped ball B, and the move involves the following two tetrahedra:
    --> the tetrahedron T that forms B; and
    --> the tetrahedron that is glued to T along the triangle opposite the
        vertex number edgeEnd of the given edge.

    This routine directly modifies the triangulation that contains the given
    edge. If an edge is currently numbered i in the triangulation, then it
    will be numbered r[i] in the triangulation after the requested 2-1 move
    has been performed (provided the move is actually legal); also, we will
    have r[i] == -1, where i is the index of the given edge (since this edge
    gets removed), and we will have r[-1] == j, where j is the index of the
    newly-created edge. If the move is not legal, the triangulation is left
    untouched and this routine returns None.

    If the triangulation containing the given edge is currently oriented,
    then this orientation will be preserved by the requested 2-1 move.
    """
    tri = edge.triangulation()

    # Is the requested 2-1 move legal?
    if not tri.twoOneMove( edge, edgeEnd, True, False ):
        return None

    # Perform the 2-1 move as a 2-3 move followed by a 2-0 move.
    emb = edge.embedding(0)
    f23 = emb.tetrahedron().triangle( emb.vertices()[edgeEnd] )
    e20 = edge.index()
    rInc = twoThree(f23)
    rDec = twoZero( tri.edge( rInc[e20] ) )

    # Work out how the edges were renumbered. Don't forget that we have
    # removed an edge, and also created an edge.
    renum = {}
    for e in range( -1, tri.countEdges() + 1 ):
        renum[e] = rDec[ rInc[e] ]
    return renum


# Tests.
if __name__ == "__main__":
    from sys import argv, stdout
    from timeit import default_timer
    from test import parseTestNames, doTest, allTestsPassedMessage

    RandomEngine.reseedWithHardware()
    #TODO Add 2-1 move tests?
    availableTests = [ "23", "20", "44" ]
    #TODO Reinstate full suite of tests once we have finished updating all the
    #       implementations.
#    availableTests = [ "23", "20", "44", "graph" ]
    testNames = parseTestNames( argv[1:], availableTests )

    #NOTE These tests use the Isomorphism3 bracket operator, introduced in
    #       Regina 7.1, to apply an isomorphism to a Triangulation3 object. In
    #       older versions of Regina, equivalent functionality was provided by
    #       the Isomorphism3.apply() routine.

    def generateIsomorphisms( triSize, maxIsos ):
        iso = Isomorphism3.identity(triSize)
        iso.inc()
        count = 0
        while ( count < maxIsos ) and ( not iso.isIdentity() ):
            count += 1
            yield iso
            for _ in range(count*25):
                iso.inc()
        return

    def multiplicities(edge):
        """
        Returns a list of the multiplicities of the given edge, sorted from
        largest to smallest.

        The multiplicity of an edge e with respect to a tetrahedron T is the
        number of times (max 6) in which e is embedded as an edge of T. The
        returned list will contain all nonzero multiplicities of the given
        edge.

        Note that the sum of the multiplicities is precisely the degree.
        """
        tri = edge.triangulation()
        ans = []
        for tet in tri.tetrahedra():
            mult = 0
            for e in range(6):
                if tet.edge(e) == edge:
                    mult += 1
            if mult:
                ans.append(mult)
        ans.sort( reverse=True )
        return ans

    def relativeEdgeOrientations( tri, iso ):
        """
        Returns the relative edge orientations resulting from applying the given
        isomorphism to the given triangulation.
        """
        newTri = iso(tri)
        ans = []
        for e in tri.edges():
            tetIndex = e.front().tetrahedron().index()
            oldVerts = e.front().vertices()
            tetImage = iso.simpImage(tetIndex)
            newVerts = iso.facetPerm(tetIndex) * oldVerts

            # Compare newVerts with the vertex permutation given by edgeMapping().
            compareVerts = newTri.tetrahedron(tetImage).edgeMapping(
                    Edge3.faceNumber(newVerts) )
            if compareVerts[0] == newVerts[0]:
                # Same orientation.
                ans.append(1)
            else:
                # Opposite orientation.
                ans.append(-1)
        return ans

    def verifyRenum( before, renum, inter, innum, after ):
        """
        Checks that renum and innum appropriately track edge renumberings
        resulting respectively from a move taking before to inter, and
        then a supposed inverse move taking inter to after.
        """
        # The idea is to search for the isomorphism that suitably relates
        # the before and after triangulations.
        moveRelOr = []
        for i in range( before.countEdges() ):
            # Relative orientation after initial move.
            tetIndex = renum[i].tetrahedron().index()
            verts = renum[i].vertices()
            compareVerts = inter.tetrahedron(tetIndex).edgeMapping(
                    Edge3.faceNumber(verts) )
            if verts[0] == compareVerts[0]:
                moveRelOr.append(1)
            else:
                moveRelOr.append(-1)

            # Relative orientation after inverse move.
            ii = edgeIndFromEmb( renum[i] )
            tetIndex = innum[ii].tetrahedron().index()
            verts = innum[ii].vertices()
            compareVerts = after.tetrahedron(tetIndex).edgeMapping(
                    Edge3.faceNumber(verts) )
            if verts[0] != compareVerts[0]:
                moveRelOr[-1] *= -1
        for iso in before.findAllIsomorphisms(after):
            if moveRelOr == relativeEdgeOrientations( before, iso ):
                # Found the desired isomorphism.
                return True
        return False

    # Test 2-3 and 3-2 moves.
    if "23" in testNames:
        print( "+-------------------+" )
        print( "| 2-3 and 3-2 moves |")
        print( "+-------------------+" )

        def test23single( face, expectedTri ):
            """
            Perform a 2-3 move on the given face, and check (among other
            things) that the result is isomorphic to expectedTri.
            """
            origTri = face.triangulation()
            tri23 = Triangulation3(origTri)
            renum = twoThree( tri23.triangle( face.index() ) )
            if not expectedTri.isIsomorphicTo(tri23):
                print(expectedTri)
                print(tri23)
                msg = "Face {}: Not isomorphic!"
                raise AssertionError( msg.format( face.index() ) )
            if ( origTri.isOriented() ) and ( not tri23.isOriented() ):
                print(tri23)
                msg = "Face {}: Failed to preserve orientation!"
                raise AssertionError( msg.format( face.index() ) )

            # Also make sure that the inverse 3-2 move gives the correct
            # triangulation, up to isomorphism.
            inv = Triangulation3(tri23)
            removedEdgeIndex = edgeIndFromEmb( renum[-1] )
            innum = threeTwo( inv.edge(removedEdgeIndex) )
            if not origTri.isIsomorphicTo(inv):
                # This test is subsumed by the more detailed isomorphisms
                # tests below, but we keep it anyway.
                print(inv)
                print( inv.detail() )
                msg = "{}: Inverse not isomorphic!"
                raise AssertionError( msg.format( face.index() ) )
            if ( tri23.isOriented() ) and ( not inv.isOriented() ):
                print(inv)
                msg = "Face {}: Inverse failed to preserve orientation!"
                raise AssertionError( msg.format( face.index() ) )
            if innum[removedEdgeIndex] is not None:
                msg = "Face {}: Inverse continues to track removed edge!"
                raise AssertionError( msg.format( face.index() ) )

            # Check that the renumberings are sensible by comparing edge
            # multiplicities in origTri and inv.
            #
            # This test is subsumed by the more detailed isomorphisms tests
            # below, but we keep it anyway.
            for i in range( origTri.countEdges() ):
                mults = multiplicities( origTri.edge(i) )
                comMults = multiplicities(
                        inv.edge(
                            edgeIndFromEmb( innum[
                                edgeIndFromEmb( renum[i] ) ] ) ) )
                if mults != comMults:
                    print( { k: edgeIndFromEmb(v)
                            for k, v in renum.items() } )
                    print( { k: edgeIndFromEmb(v)
                            for k, v in innum.items() if v is not None } )
                    msg = "Face {}: Unmatched edge multiplicities!"
                    raise AssertionError( msg.format( face.index() ) )

            # Check that the renumberings are sensible.
            if not verifyRenum( origTri, renum, tri23, innum, inv ):
                print( { k: edgeIndFromEmb(v)
                        for k, v in renum.items() } )
                print( { k: edgeIndFromEmb(v)
                        for k, v in innum.items() if v is not None } )
                msg = "Face {}: Relabellings failed!"
                raise AssertionError( msg.format( face.index() ) )

            # All done!
            return

        def test23all( testSig, maxIsos=16 ):
            print( "2-3 and 3-2 moves on \"{}\"".format(testSig) )
            stdout.flush()
            t = Triangulation3.fromIsoSig(testSig)
            t.orient()

            # Test 2-3 moves on all eligible triangles.
            count = 0
            for f in t.triangles():
                # Make sure that we get the correct isomorphism type, and that
                # we preserve orientation.
                #
                #NOTE Triangulation3.withPachner(f) was introduced in
                #       Regina 7.4. Older versions of Regina did not provide
                #       a routine to construct a new triangulation via a 2-3
                #       move; instead, this behaviour was achieved by first
                #       constructing a copy of the triangulation, and then
                #       performing the appropriate 2-3 move one the copy.
                pach = t.withPachner(f)
                if pach is None:
                    # 2-3 move not eligible on f.
                    continue
                count += 1
                try:
                    test23single( f, pach )
                except AssertionError as ae:
                    print(t)
                    raise ae

                # To test as many cases of the implementation as possible,
                # test the same 2-3 move with several relabellings of t.
                for iso in generateIsomorphisms( t.size(), maxIsos ):
                    r = iso(t)
                    source = f.embedding(0).tetrahedron().index()
                    fnum = f.embedding(0).face()
                    tetImage = iso.simpImage(source)
                    faceImage = iso.facetPerm(source)[fnum]
                    try:
                        test23single(
                                r.tetrahedron(tetImage).triangle(faceImage),
                                pach )
                    except AssertionError as ae:
                        print(iso)
                        raise ae

            # All done!
            print( "Tested {} pairs of 2-3 and 3-2 moves.".format(count) )
            print()
            stdout.flush()
            return

        # Run 2-3 and 3-2 move tests.
        start = default_timer()
        for testSig in [ "cMcabbgqs", "gLLPQcdefeffpvauppb" ]:
            test23all(testSig)
        print()
        print( "Time: {:.6f}".format( default_timer() - start ) )
        print( "2-3 and 3-2 moves: All tests passed!" )
        print()
        stdout.flush()

    # Test 2-0 moves.
    if "20" in testNames:
        print( "+-----------+" )
        print( "| 2-0 moves |")
        print( "+-----------+" )

        def zeroTwo( tri, edgeIndex, i, ii ):
            """
            Performs a 0-2 move and returns the corresponding edge-renumbering
            map, or None if the requested 0-2 move is not legal.
            """
            #NOTE This implementation assumes that the two new tetrahedra
            #       introduced by move02() are located at the last 2 indices.

            # Tracking edge embeddings through the requested 0-2 move is easy
            # because every existing tetrahedron survives the move.
            relab = EdgeLabelling(tri)
            if not tri.move02( tri.edge(edgeIndex), i, ii ):
                return None

            # Find the new edge that the 0-2 move just created.
            tet = tri.tetrahedron( tri.size() - 1 )
            for edgeNum in range(6):
                newEdge = tet.edge(edgeNum)
                if newEdge.degree() != 2:
                    continue

                # Because the last two tetrahedra of tri were introduced by
                # the 0-2 move that we just performed, the newly-introduced
                # edge must be the unique degree-2 edge that is incident to
                # both of the last two tetrahedra (this is easy to check).
                incidentTetInds = { newEdge.front().tetrahedron().index(),
                                   newEdge.back().tetrahedron().index() }
                if incidentTetInds == { tri.size() - 1, tri.size() - 2 }:
                    relab[-1] = newEdge.front()
                    return relab
            raise AssertionError(
                    "zeroTwo() should never reach this point." )

        def test20single( tri, edgeIndex, i, ii ):
            """
            Checks that the inverse 2-0 move correctly inverts the given 0-2
            move.
            """
            tri02 = Triangulation3(tri)
            renum = zeroTwo( tri02, edgeIndex, i, ii )
            if renum is None:
                # No 0-2 move here.
                return False
            msg = "Inverse 2-0 move {}, {}, {}: ".format( edgeIndex, i, ii )

            # Test that the inverse 2-0 move, performed on tri02 using
            # twoZero(), brings us back to a triangulation isomorphic to tri.
            inv = Triangulation3(tri02)
            innum = twoZero( inv.edge( edgeIndFromEmb( renum[-1] ) ) )
            if not inv.isIsomorphicTo(tri):
                print(tri)
                print(tri02)
                print(inv)
                raise AssertionError( msg + "Not isomorphic!" )
            if ( tri02.isOriented() ) and ( not inv.isOriented() ):
                print(tri02)
                print(inv)
                raise AssertionError( msg + "Failed to preserve orientation!" )

            # Check that the renumberings are sensible.
            if not verifyRenum( tri, renum, tri02, innum, inv ):
                print( { k: edgeIndFromEmb(v)
                        for k, v in renum.items() } )
                print( { k: edgeIndFromEmb(v)
                        for k, v in innum.items() if v is not None } )
                raise AssertionError( msg + "Relabellings failed!" )

            # All done!
            return True

        def test20all(testSig):
            """
            Test inverse 2-0 moves on the triangulations obtained by
            performing 0-2 moves on the given iso sig.
            """
            print( "0-2 and 2-0 moves on \"{}\"".format(testSig) )
            stdout.flush()
            t = Triangulation3.fromIsoSig(testSig)
            t.orient()

            count = 0
            for e in t.edges():
                deg = e.degree()
                for i in range(deg):
                    for ii in range( i, deg ):
                        if test20single( t, e.index(), i, ii ):
                            count += 1

            # All done!
            print( "Tested {} 2-0 moves.".format(count) )
            print()
            stdout.flush()
            return

        # Run 2-0 move tests.
        start = default_timer()
        test20all( "gLLPQcdefeffpvauppb" )
        print()
        print( "Time: {:.6f}".format( default_timer() - start ) )
        print( "2-0 moves: All tests passed!" )
        print()
        stdout.flush()

    # Test 4-4 moves.
    if "44" in testNames:
        print( "+-----------+" )
        print( "| 4-4 moves |")
        print( "+-----------+" )

        def test44single(edge):
            """
            Tests all possible 4-4 moves on the given edge.
            """
            t = edge.triangulation()
            #NOTE Triangulation3.has44( e, ax ) was introduced in Regina 7.4.
            #       In older versions of Regina, equivalent functionality
            #       (checking eligibility of the move, but not performing it)
            #       was provided by
            #           Triangulation3.fourFourMove( e, ax, True, False ).
            if not t.has44( edge, 0 ):
                # A 4-4 move with newAxis == 0 is available if and only if a
                # 4-4 move with newAxis == 1 is available.
                return False

            # The input edge forms one of three possible axes for an
            # octahedron built from four tetrahedra. At each such axis, we
            # perform both possible 4-4 moves, and check that the isomorphism
            # types of the resulting triangulations all match up.
            isoSigSet = { i: set() for i in range(3) }
            isoSigSet[2].add( t.isoSig() )
            for newAxis in range(2):
                #NOTE Triangulation3.with44( e, ax ) was introduced in
                #       Regina 7.4. Older versions of Regina did not provide a
                #       routine to construct a new triangulation via a 4-4
                #       move; instead, this behaviour was achieved by first
                #       constructing a copy of the triangulation, and then
                #       performing the appropriate 4-4 move on the copy.
                reg44 = t.with44( edge, newAxis )
                new44 = Triangulation3(t)
                relab = fourFour( new44.edge( edge.index() ), newAxis )

                # Test that fourFour gives the right isomorphism type.
                move = "4-4 move {}/{}: ".format( edge.index(), newAxis )
                if not reg44.isIsomorphicTo(new44):
                    print(t)
                    print(reg44)
                    print(new44)
                    raise AssertionError( move + " Not isomorphic!" )
                newSig = new44.isoSig()
                isoSigSet[2].add(newSig)
                isoSigSet[newAxis].add(newSig)

                # Sanity checks on the relabelling.
                if relab[ edge.index() ] is not None:
                    raise AssertionError(
                            move +
                            " Relabelling continues tracking removed edge!" )
                if relab[-1] is None:
                    raise AssertionError(
                            move +
                            " Relabelling fails to track new edge!" )

                # Finish filling in isoSigSet so that we can test the 4-4
                # moves on the new edge that we just created.
                for invAxis in range(2):
                    inv44 = Triangulation3(new44)
                    fourFour( inv44.edge( edgeIndFromEmb( relab[-1] ) ),
                             invAxis )
                    isoSigSet[newAxis].add( inv44.isoSig() )
            for newAxis in range(2):
                if isoSigSet[2] != isoSigSet[newAxis]:
                    raise AssertionError(
                            "Edge {}: Error with newAxis {}.".format(
                                edge.index(), newAxis ) )
            return True

        def test44all(testSig):
            print( "4-4 moves on \"{}\"".format(testSig) )
            stdout.flush()
            t = Triangulation3.fromIsoSig(testSig)
            t.orient()

            count = 0
            for e in t.edges():
                if test44single(e):
                    count += 1
            print( "Tested 4-4 moves on {} edges.".format(count) )
            print()
            stdout.flush()
            return

        # Run 4-4 move tests.
        start = default_timer()
        for testSig in [ "gLLPQcdefeffpvauppb",
                        "gLLPQceeffefiiaealx",
                        "gvLQQcdeffeffffaafa",
                        "gLLAQcdcdfffpvbbbvo",
                        "ivLAPQcdefeghghhbbpbuabbv" ]:
            test44all(testSig)
        print()
        print( "Time: {:.6f}".format( default_timer() - start ) )
        print( "4-4 moves: All tests passed!" )
        print()
        stdout.flush()

    #TODO Update tests to use new EdgeLabelling.

    # Perform a bunch of moves to check that we don't get any exceptions.
    if "graph" in testNames:
        print( "+---------------------+" )
        print( "| Triangulation graph |" )
        print( "+---------------------+" )

        print()
        print( "Perform moves on \"cMcabbgqs\" up to height 3." )
        print()
        initSig = "cMcabbgqs"
        stack = [ initSig ]
        seen = { initSig }
        initSize = 2
        maxSize = 5
        counts = { "2-1": 0, "2-0": 0, "3-2": 0, "4-4": 0, "2-3": 0 }
        while stack:
            sig = stack.pop()
            tri = Triangulation3.fromIsoSig(sig)
            if tri.size() <= initSize:
                print(sig)
                stdout.flush()

            # Moves on edges.
            for e in tri.edges():
                if e.degree() == 1:
                    for edgeEnd in {0,1}:
                        newTri = Triangulation3(tri)
                        renum = twoOne( newTri.edge( e.index() ), edgeEnd )
                        if renum is None:
                            continue
                        counts["2-1"] += 1

                        # If we haven't seen the new triangulation before (up to
                        # combinatorial isomorphism), then add it to the stack.
                        newSig = newTri.isoSig()
                        if newSig not in seen:
                            seen.add(newSig)
                            stack.append(newSig)
                elif e.degree() == 2:
                    newTri = Triangulation3(tri)
                    renum = twoZero( newTri.edge( e.index() ) )
                    if renum is None:
                        continue
                    counts["2-0"] += 1

                    # If we haven't seen the new triangulation before (up to
                    # combinatorial isomorphism), then add it to the stack.
                    newSig = newTri.isoSig()
                    if newSig not in seen:
                        seen.add(newSig)
                        stack.append(newSig)
                elif e.degree() == 3:
                    newTri = Triangulation3(tri)
                    renum = threeTwo( newTri.edge( e.index() ) )
                    if renum is None:
                        continue
                    counts["3-2"] += 1

                    # If we haven't seen the new triangulation before (up to
                    # combinatorial isomorphism), then add it to the stack.
                    newSig = newTri.isoSig()
                    if newSig not in seen:
                        seen.add(newSig)
                        stack.append(newSig)
                elif e.degree() == 4:
                    for newAxis in {0,1}:
                        newTri = Triangulation3(tri)
                        renum = fourFour( newTri.edge( e.index() ), newAxis )
                        if renum is None:
                            continue
                        counts["4-4"] += 1

                        # If we haven't seen the new triangulation before (up to
                        # combinatorial isomorphism), then add it to the stack.
                        newSig = newTri.isoSig()
                        if newSig not in seen:
                            seen.add(newSig)
                            stack.append(newSig)

            # 2-3 moves (only if such moves would not exceed the given height).
            if tri.size() == maxSize:
                continue
            for f in tri.triangles():
                newTri = Triangulation3(tri)
                renum = twoThree( newTri.triangle( f.index() ) )
                if renum is None:
                    continue
                counts["2-3"] += 1

                # If we haven't seen the new triangulation before (up to
                # combinatorial isomorphism), then add it to the stack.
                newSig = newTri.isoSig()
                if newSig not in seen:
                    seen.add(newSig)
                    stack.append(newSig)

        # Print move counts.
        print()
        print( "Tested: {} {}; {} {}; {} {}; {} {}; {} {}.".format(
            counts["2-1"], "2-1 moves",
            counts["2-0"], "2-0 moves",
            counts["3-2"], "3-2 moves",
            counts["4-4"], "4-4 moves",
            counts["2-3"], "2-3 moves" ) )
