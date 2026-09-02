"""
Scratch work for experimental testing of recogniseSFS().

Usage: Must supply the following integer input arguments:
    (1) Number of workers
    (2) Signed genus of the base surface
    (3) Number of boundary components
    (4) Number of exceptional fibres
    (5) Maximum multiplicity of fibres
    (6) Maximum number of randomisation attempts for each SFS
"""
import sys
from math import gcd as pythonGCD
from itertools import combinations_with_replacement as combinsWithRep
from multiprocessing import Pool, TimeoutError
from timeit import default_timer
from regina import *
from construct.sfs import orientableSFS


def attemptHardSFS( baseSignedGenus, boundaries, fibres, attempts ):
    """
    Attempts to construct a hard triangulation of an orientable Seifert fibre
    space by randomising a few times.

    If one such randomisation attempt succeeds, then this routine returns a
    2nd-generation isomorphism signature for the hard triangulation.

    Here, we consider a triangulation "hard" if, after simplifying using
    Triangulation3.simplify(), combinatorial recognition via Regina's
    BlockedSFS.recognise() function fails.

    In detail, the arguments to this routine should specify such a Seifert
    fibre space as follows:
    --> baseSignedGenus should be an integer giving the "signed genus" of the
        base surface. That is:
        --- if baseSignedGenus >= 0, then the base will be orientable and
            have genus baseSignedGenus; and
        --- if baseSignedGenus < 0, then the base will be the surface with
            nonorientable genus -baseSignedGenus.
    --> boundaries should be a non-negative integer specifying the number of
        boundary components of the Seifert fibre space.
    --> fibres should be a collection of exceptional fibres, where each fibre
        is specified by a pair (a,b) of integers such that:
        --- a >= 1; and
        --- a and b are coprime.

    The attempts parameter specifies the maximum number of randomisation
    attempts. If this is negative, then this routine will attempt an
    unbounded number of randomisations (and hence might *never* terminate).
    Otherwise, if this is non-negative, and if none of the attempts leads to
    a hard triangulation, then this routine will simply return None.
    """
    tri = orientableSFS( baseSignedGenus, boundaries, *fibres )
    max23 = 512
    alwaysModify = True
    i = 0
    while i != attempts:
        i += 1
        try:
            tri.simplifyUpDown( max23, alwaysModify )
        except AttributeError:
            raise RuntimeError( "attemptHardSFS() requires a sufficiently " +
                               "new version of Regina that includes the " +
                               "Triangulation3.simplifyUpDown() routine" )

        # Try really hard to simplify.
        tri.simplify()
        simplifiedNow = True    # Attempt simplification at least twice
        while simplifiedNow:
            simplifiedNow = tri.simplify()

        # Is the randomised triangulation hard?
        blocked = BlockedSFS.recognise(tri)
        if blocked is None:
            return tri.neoSig()
    return None


if __name__ == "__main__":
    numWorkers = int( sys.argv[1] )
    baseSignedGenus = int( sys.argv[2] )
    numBdries = int( sys.argv[3] )
    if numBdries < 1:
        raise ValueError( "Can't have {} boundaries. ".format(numBdries) +
                         "Must be at least 1." )
    numFibres = int( sys.argv[4] )
    if numFibres < 0:
        raise ValueError( "Can't have {} ".format(numFibres) +
                         "exceptional fibres. Must be at least 0." )
    maxMultiplicity = int( sys.argv[5] )
    if maxMultiplicity < 2:
        raise ValueError( "Can't have exceptional fibres of " +
                         "multiplicity {}. ".format(maxMultiplicity) +
                         "Must be at least 2." )
    maxAttempts = int( sys.argv[6] )
    if baseSignedGenus >= 0:
        orbl = "Orientable"
        baseClass = SFSpace.Class.bo1
        baseGenus = baseSignedGenus
    else:
        orbl = "Non-orientable"
        baseClass = SFSpace.Class.bn2
        baseGenus = -baseSignedGenus
    print( "Workers: {}".format(numWorkers) )
    print( "{} base, genus {} with {} boundaries.".format(
        orbl, baseGenus, numBdries ) )
    print( "{} fibres with max multiplicity {}.".format(
        numFibres, maxMultiplicity ) )
    if maxAttempts < 0:
        print( "Unlimited number of randomisation attempts." )
    else:
        print( "Max {} randomisation attempts per SFS.".format(maxAttempts) )
    print()
    sys.stdout.flush()
    start = default_timer()
    fibres = []
    for p in range( 2, 1 + maxMultiplicity ):
        for q in range( 1, p ):
            if pythonGCD(p, q) != 1:
                continue
            fibres.append( (p, q) )

    # Try to generate a non-standard triangulation for each Seifert fibred
    # space satisfying the given parameters.
    unfinished = dict()
    manifold = dict()
    with Pool( processes=numWorkers ) as pool:
        for i, fibreParams in enumerate(
                combinsWithRep( fibres, numFibres ) ):
            #TODO Filter out fibre parameters that are just negations
            unfinished[i] = pool.apply_async( attemptHardSFS, args = (
                baseSignedGenus, numBdries, fibreParams, maxAttempts ) )
            manifold[i] = SFSpace( baseClass, baseGenus, numBdries )
            for p, q in fibreParams:
                manifold[i].insertFibre( p, q )
        print( "Generating up to {} triangulations".format(i + 1) )
        print()
        sys.stdout.flush()

        # Repeatedly poll for new results.
        found = 0
        while unfinished:
            nowFinished = set()
            for i in unfinished:
                try:
                    output = unfinished[i].get( timeout=0.001 )
                except TimeoutError:
                    pass
                else:
                    # Got a new result!
                    nowFinished.add(i)

                    # Did we actually find a hard SFS triangulation?
                    if output is not None:
                        found += 1
                        print( "#{}. Time: {:.6f}. {}. \"{}\".".format(
                            found, default_timer() - start, manifold[i], output ) )
            for i in nowFinished:
                del unfinished[i]
