"""
Blueprints for reconstructing 3-manifold triangulations.

These blueprints serve two functions:
--> They are picklable, which facilitates passing the data of a triangulation
    between processes (at the cost of having to first compute the blueprint,
    and later reconstructing the triangulation).
--> Unlike isomorphism signatures, which only reconstruct triangulations up to
    combinatorial isomorphism, these blueprints keep track of all the
    tetrahedron indices and face numberings.
"""
from regina import *


def triangulationBlueprint(tri):
    """
    Returns a picklable blueprint for the given 3-manifold triangulation tri.

    In detail, this routine returns a tuple (S,T,F), where:
    --> S is the isomorphism signature for tri;
    --> T is a list such that T[i] gives the index of the tetrahedron in
        newTri := Triangulation3.fromIsoSig(S) that corresponds to
        tetrahedron i of tri; and
    --> F is a list such that F[i] is length-4 list with F[i][j] being
        the vertex number of newTri.tetrahedron( T[i] ) that corresponds
        to vertex j of tri.tetrahedron(i).
    The returned blueprint can be used, via the reconstructTriangulation()
    routine, to build a clone of the given triangulation. Specifically, the
    following command would return a clone of tri:
        reconstructTriangulation( *triangulationBlueprint(tri) )
    """
    sig, isom = tri.isoSigDetail()

    # Convert the isomorphism into the lists T and F.
    tetImages = []
    facePerms = []
    for i in range( tri.size() ):
        tetImages.append( isom.simpImage(i) )
        facePerms.append(
                [ isom.facetPerm(i)[j] for j in range(4) ] )
    return ( sig, tetImages, facePerms )


def reconstructTriangulation( sig, tetImages, facePerms ):
    """
    Returns a triangulation reconstructed from the given data.

    The main purpose of this routine is to construct a clone of a 3-manifold
    triangulation from a collection of picklable data. As described in the
    documentation for the triangulationBlueprint() routine, the following
    command would return a clone of a 3-manifold triangulation tri:
        reconstructTriangulation( *triangulationBlueprint(tri) )
    """
    tri = Triangulation3.fromIsoSig(sig)

    # Adjust tri to use the original labelling (instead of the canonical
    # labelling given by fromIsoSig()).
    isom = Isomorphism3( tri.size() )
    for i in range( tri.size() ):
        isom.setSimpImage( i, tetImages[i] )
        isom.setFacePerm( i, Perm4( *facePerms[i] ) )
    tri = isom.inverse()(tri)
    return tri
