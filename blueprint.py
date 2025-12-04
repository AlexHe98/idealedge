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

    In detail, this routine returns a pair (S,G), where:
    --> S is the size (i.e., the number of tetrahedra) of tri.
    --> G is a list such that each entry is of the form [i,f,j,p], and
        describes a gluing of tri as follows:
        --- i is a tetrahedron index of tri;
        --- f is a face number of tetrahedron i;
        --- j is the index of the tetrahedron adjacent to tetrahedron i along
            face f; and
        --- p specifies that the gluing along face f of tetrahedron i is given
            by the permutation Perm4.S4[p].
    This is essentially the data required to construct a triangulation using
    the Triangulation3.fromGluings() routine, except that we use integer
    indices (denoted p above) to specify the gluing permutations using a
    picklable data type.

    The returned blueprint can be used, via the reconstructionTriangulation()
    routine, to build a clone of the given triangulation. Specifically, the
    following command would return a clone of tri:
        reconstructTriangulation( *triangulationBlueprint(tri) )
    """
    return ( tri.size(), _gluings(tri) )


def reconstructTriangulation( size, gluings ):
    """
    Returns a triangulation reconstructed from the given data.

    The main purpose of this routine is to construct a clone of a 3-manifold
    triangulation from a collection of picklable data. As described in the
    documentation for the triangulationBlueprint() routine, the following
    command would return a clone of a 3-manifold triangulation tri:
        reconstructTriangulation( *triangulationBlueprint(tri) )
    """
    # Convert the picklable gluings into a description that we can pass to
    # Triangulation3.fromGluings().
    use = []
    for g in gluings:
        use.append( [ g[0], g[1], g[2], Perm4.S4[g[3]] ] )
    return Triangulation3.fromGluings( size, use )


def _gluings(tri):
    """
    Returns a picklable description of the gluings used to construct tri.
    """
    # We could get this information by parsing the output of tri.source(),
    # but it seems better to just compute this for ourselves.
    ans = []
    for face in tri.triangles():
        # By iterating through the triangles, we ensure that we only encode
        # each gluing once.
        myTet = face.front().tetrahedron()
        myFace = face.front().face()
        yourTet = myTet.adjacentTetrahedron(myFace)
        if yourTet is None:
            continue
        gluing = myTet.adjacentGluing(myFace)
        ans.append( [ myTet.index(), myFace, yourTet.index(),
                     gluing.S4Index() ] )
    return ans
