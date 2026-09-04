"""
Scratch work for experimental testing of recogniseSFS().

Usage: Must supply the following integer input arguments:
    (1) Signed genus of the base surface
    (2) Number of boundary components
    (3) Number of exceptional fibres
    (4) Maximum multiplicity of fibres
    (5) Maximum number of randomisation attempts for each SFS
"""
import sys
from math import gcd as pythonGCD
from itertools import combinations_with_replacement as combinsWithRep
from multiprocessing import Pool, TimeoutError, cpu_count
from timeit import default_timer
from collections import Counter
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


def _hardSFSPath( isBaseOrbl, genus, numBdries ):
    if isBaseOrbl:
        base = "baseOr"
    else:
        base = "baseNor"
    return "sfs-samples/sfs-{}-g{}-m{}.txt".format( base, genus, numBdries )


def _sortedFibreParams(fibreParams):
    return tuple( sorted(fibreParams) )


def _sortedFibreNegation(fibreParams):
    return _sortedFibreParams(
            (2, 1) if p == 2 else (p, -q) for p, q in fibreParams )


def _readSFSLine(line):
    data = line.rstrip().split(" ")
    fibres = []
    neoSig = data[0]
    for fibreString in data[1:]:
        fibres.append( tuple(
            int(k) for k in fibreString[1:-1].split(",") ) )
    return neoSig, tuple(fibres)


def _sfsLine( neoSig, fibres ):
    line = neoSig
    for p, q in fibres:
        line += " ({},{})".format( p, q )
    return line + "\n"


if __name__ == "__main__":
    baseSignedGenus = int( sys.argv[1] )
    numBdries = int( sys.argv[2] )
    if numBdries < 1:
        raise ValueError( "Can't have {} boundaries. ".format(numBdries) +
                         "Must be at least 1." )
    numFibres = int( sys.argv[3] )
    if numFibres < 0:
        raise ValueError( "Can't have {} ".format(numFibres) +
                         "exceptional fibres. Must be at least 0." )
    maxMultiplicity = int( sys.argv[4] )
    if maxMultiplicity < 2:
        raise ValueError( "Can't have exceptional fibres of " +
                         "multiplicity {}. ".format(maxMultiplicity) +
                         "Must be at least 2." )
    maxAttempts = int( sys.argv[5] )
    isBaseOrbl = ( baseSignedGenus >=  0 )
    if isBaseOrbl:
        orbl = "Orientable"
        baseClass = SFSpace.Class.bo1
        baseGenus = baseSignedGenus
    else:
        orbl = "Non-orientable"
        baseClass = SFSpace.Class.bn2
        baseGenus = -baseSignedGenus
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
    filePath = _hardSFSPath( isBaseOrbl, baseGenus, numBdries )
    with open( filePath, 'a+' ) as sfsFile:
        # First read in all the manifolds where we have previously found
        # non-standard triangulations.
        foundParams = set()
        sortedNegations = set()
        sfsFile.seek(0)
        for line in sfsFile:
            _, fibres = _readSFSLine(line)
            foundParams.add(fibres)
            sortedNegations.add( _sortedFibreNegation(fibres) )

        # For each SFS satisfying the given parameters, if we haven't already
        # found a non-standard triangulation, then try to generate one now.
        start = default_timer()
        fibres = []
        for p in range( 2, 1 + maxMultiplicity ):
            for q in range( 1, 1 + p//2 ):
                if pythonGCD(p, q) != 1:
                    continue
                fibres.append( (p, q) )
                if p != 2:
                    fibres.append( (p, -q) )
        paramsToAttempt = dict()
        manifold = dict()
        count = 0
        for fibreParams in combinsWithRep( fibres, numFibres ):
            if fibreParams in foundParams:
                continue

            # Avoid duplicates with the opposite orientation.
            if _sortedFibreParams(fibreParams) in sortedNegations:
                continue
            sortedNegations.add( _sortedFibreNegation(fibreParams) )

            # We have a bounded orientable SFS for which we hope to find
            # a non-standard triangulation.
            paramsToAttempt[count] = fibreParams
            manifold[count] = SFSpace( baseClass, baseGenus, numBdries )
            for p, q in fibreParams:
                manifold[count].insertFibre( p, q )
            count += 1

        # Do as many in parallel as possible.
        numWorkers = min( count, max( 1, cpu_count() - 1 ) )
        print( "Time: {:.6f}".format( default_timer() - start ) )
        print( "Generating up to {} new triangulations".format(count) )
        print( "Using {} workers".format(numWorkers) )
        print()
        sys.stdout.flush()
        unfinished = dict()
        with Pool( processes=numWorkers ) as pool:
            for i in range(count):
                unfinished[i] = pool.apply_async( attemptHardSFS, args = (
                    baseSignedGenus, numBdries,
                    paramsToAttempt[i], maxAttempts ) )

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
                            # Write to file and print update to stdout.
                            found += 1
                            sfsFile.write( _sfsLine(
                                output, paramsToAttempt[i] ) )
                            sfsFile.flush()
                            print( "#{}. Time: {:.6f}. {}. \"{}\".".format(
                                found, default_timer() - start, manifold[i], output ) )
                            sys.stdout.flush()
                for i in nowFinished:
                    del unfinished[i]
    print()
    print( "Time: {:.6f}".format( default_timer() - start ) )
