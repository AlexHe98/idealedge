"""
Scratch work for bounded orientable SFS recognition using edge-ideal
triangulations.
"""
import sys
from timeit import default_timer
from regina import *
from construct.sfs import orientableSFS
from recsfs import recogniseSFS


#TODO Make this more general.
def filledHomology(annulus):
    """
    Given an annulus A in a triangulation whose boundary is a two-triangle
    torus B, such that both boundary curves of A are parallel to an edge of B,
    computes the first homology of the manifold given by filling along the
    slope of the boundary of A.
    """
    # What's the homology of the manifold before filling?
    tri = annulus.triangulation()
    markedH1 = tri.markedHomology(1)
    rank = markedH1.rank()
    numInvFacs = markedH1.countInvariantFactors()
    invFacs = [ markedH1.invariantFactor(i) for i in range(numInvFacs) ]
    numGenerators = rank + numInvFacs
    # Each column represents a generator.
    # Each row gives a relation.
    presentation = []
    for i in range(numInvFacs):
        relation = [0] * numGenerators
        relation[i] = invFacs[i]
        presentation.append(relation)

    # We need to add one extra relation corresponding to the curve c that is
    # killed by the Dehn filling. By assumption, there is a boundary edge that
    # is parallel to c.
    missedEdgeIndex = None
    for e in tri.edges():
        if not e.isBoundary():
            continue
        if annulus.edgeWeight( e.index() ).pythonValue() == 0:
            missedEdgeIndex = e.index()
            break
    if missedEdgeIndex is None:
        raise ValueError( "Given surface doesn't satisfy preconditions" )

    # Write the missed edge in terms of generators, so that we can then add a
    # relation that kills c.
    cycle = [0] * tri.countEdges()
    cycle[missedEdgeIndex] = 1
    relation = markedH1.snfRep(cycle)
    #TODO BEGIN TEST
    if presentation:
        print( AbelianGroup( MatrixInt(presentation) ) )
    print(relation)
    #TODO END TEST
    presentation.append(relation)
    return AbelianGroup( MatrixInt(presentation) )


if __name__ == "__main__":
    genus = int( sys.argv[1] )
    boundaries = int( sys.argv[2] )
    params = [ int(n) for n in sys.argv[3:] ]
    fibres = []
    while params:
        q = params.pop()
        p = params.pop()
        fibres.append( (p,q) )
    print( "g={}, b={}".format( genus, boundaries ) )
    print(fibres)
    print()
    tri = orientableSFS( genus, boundaries, *fibres )
#    hardS3 = Triangulation3("-cwc3admabhuOIgKiJlYydn+idrcreqg5Kuazew8OJu+qLwkrexkzfuEHLwwrLzMbgyGzMyKHgB0PhKmIOFesjJkQOMkkPKGQQMO6PTOkQOWAmX66RVSYlYwdUVcRUYuRV2sdV2CJo+cCqa3ZV9YJo-0tocNlp9smrfleWcjCq-qesjfKsiF0XdF0XkLKXfBKXgBCsjPCskRSY--Vpp+0KL3-sTlHIjai+tTc8rxKyobacJ-75w-hwUWhqnMaPcGD2bXbWaXA0br8QASXWXgt1Drz3sTLirMDGAintbar1ivHubIAiytCQ62CYHqybAg7KtZzo1OigAfrN1K5qyG4ccTswDfCgBxZGmXb")
#    tri.connectedSumWith(hardS3)

    # Simplify and attempt combinatorial recognition.
    print( "Combinatorial recognition" )
    print( "-------------------------" )
    print( "Initial size: {}".format( tri.size() ) )
    sys.stdout.flush()
    start = default_timer()
    simplifiedNow = tri.simplify()
    if not simplifiedNow:
#        tri.simplifyUpDown(1024)
        # Make sure to attempt simplification at least once more.
        simplifiedNow = True
    while simplifiedNow:
        simplifiedNow = tri.simplify()
#        if not simplifiedNow:
#            simplifiedNow = tri.simplifyUpDown(1024)
    print( "Simplified size: {}".format( tri.size() ) )
    print( "Time: {:.6f}".format( default_timer() - start ) )
    blocked = BlockedSFS.recognise(tri)
    if blocked is None:
        print( "Not recognised" )
    else:
        print( "Recognised!" )
        print( blocked.manifold() )
    print( "Time: {:.6f}".format( default_timer() - start ) )
    print()
    sys.stdout.flush()

    # Now run the normal surface computation.
    print( "Normal surfaces" )
    print( "---------------" )
    sys.stdout.flush()
    start = default_timer()
    print( recogniseSFS( tri, False ) )
    print( "Time: {:.6f}".format( default_timer() - start ) )
    print()
    sys.stdout.flush()
