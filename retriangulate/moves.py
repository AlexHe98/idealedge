"""
Perform 2-3, 3-2 and 2-0 moves while tracking how edges are relabelled.
"""
from regina import *
from retriangulate.edgelabel import EdgeLabelling
from aux.tetrenum import tetRenumbering, tetIndexRenumbering


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
    newIndex = tetRenumbering(doomed)
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
    for edgeInd, emb in edgeLab.items():
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
    newIndex = tetRenumbering(doomed)
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
    for edgeInd, emb in edgeLab.items():
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
    doomedIndices = { emb.tetrahedron().index()
                     for emb in edge.embeddings() }
    newIndex = tetIndexRenumbering( tri, doomedIndices )

    # For each edge e (other than the input edge) that is tracked by the
    # reference labelling, find one tetrahedron that will meet e after we have
    # performed the requested 2-0 move.
    newEdgeLocations = dict()
    for edgeInd, emb in edgeLab.items():
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
            oldTetInd = otherEmb.tetrahedron().index()
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
    relab = twoThree( f23, edgeLab )
    return threeTwo( doomed[3].edge(e32), relab )


def twoOne( edge, edgeEnd, edgeLab=None ):
    """
    Performs a 2-1 move about the given edge, and returns an EdgeLabelling
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

    If the 2-1 move is legal, then the given edge forms the internal edge of a
    snapped ball B, and the move involves the following two tetrahedra:
    --> the tetrahedron S that forms B; and
    --> the tetrahedron that is glued to S along the triangle opposite vertex
        number edgeEnd of the given edge (hence, note that edgeEnd must be
        either 0 or 1).

    If the requested move is not legal, then edge.triangulation() is left
    entirely untouched, and this routine returns None.

    Otherwise, if the move is legal, then this routine directly modifies the
    triangulation T := edge.triangulation(). The returned EdgeLabelling r is
    structured as follows:
    --> The requested 2-1 move always destroys the given edge, so for any
        index i in the reference labelling that corresponds to the given edge,
        r[i] will be None.
    --> For every other index i in the reference labelling, which corresponds
        to some edge e in T (other than the given edge) before the move, r[i]
        will be an EdgeEmbedding3 object describing an embedding of e in T
        *after* performing the move. The embedding of e in r will have the
        same orientation as the embedding of e in the reference labelling.
    --> A 2-1 move also always creates a new edge, and r[i] will give an
        embedding of this new edge, where i is the largest negative index that
        is not already tracked by the reference labelling (typically, this
        will mean that i is -1).

    This routine might make temporary modifications to edgeLab (if supplied),
    but any such modifications will be reverted by the time this routine
    terminates.

    Note also that a 2-1 move merges two edges into a single new edge e. If
    the reference labelling tracks indices, say i and j, for both of the
    merged edges, then in the returned EdgeLabelling r, we will have that r[i]
    and r[j] give two (possibly equal, possibly distinct) EdgeEmbedding3
    objects corresponding to the same new edge e.

    If edge.triangulation() is currently oriented, then this orientation will
    be preserved by the requested 2-1 move.
    """
    tri = edge.triangulation()

    # Is the requested 2-1 move legal?
    #
    #NOTE Triangulation3.has21( e, ed ) was introduced in Regina 7.4. In older
    #       versions of Regina, equivalent functionality (checking eligibility
    #       of the move, but not performing it) was provided by
    #       Triangulation3.twoOneMove( e, ed, True, False ).
    if not tri.has21( edge, edgeEnd ):
        return None

    # To perform the move, we need to ensure that the input edge is tracked.
    if edgeLab is None:
        edgeLab = EdgeLabelling(tri)
        trackInputEdge = True
        e20 = edge.index()
    else:
        # With a custom reference labelling, we will need to work out how the
        # input edge is labelled (assuming it is even tracked).
        trackInputEdge = False
        for ei in edgeLab:
            if edge.index() != edgeLab.underlyingEdgeIndex(ei):
                continue
            trackInputEdge = True
            e20 = ei
            break
        if not trackInputEdge:
            # Temporarily track the input edge on the right (the left is
            # already reserved for tracking the newly-created edge).
            e20 = 1 + max( edgeLab.trackedIndices() )
            edgeLab[e20] = edge.front()

    # Perform the 2-1 move as a 2-3 move followed by a 2-0 move.
    emb = edge.front()
    f23 = emb.tetrahedron().triangle( emb.vertices()[edgeEnd] )
    relab = twoThree( f23, edgeLab )
    relab = twoZero(
            tri.edge( relab.underlyingEdgeIndex(e20) ),
            relab )
    if not trackInputEdge:
        relab.untrack(e20)
        # Don't forget to restore edgeLab to its original state.
        edgeLab.untrack(e20)
    return relab
