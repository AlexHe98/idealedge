"""
Auxiliary functions for working with Seifert fibrations.
"""


def sortedFibreParams(fibreParams):
    """
    Returns a sorted tuple of the given sequence of exceptional fibre
    parameters.

    Each exceptional fibre should be given as a pair (p, q) of integers such
    that p >= 2 and gcd(p, q) == 1.

    If each exceptional fibre is written in a normalised form (for example,
    for a bounded SFS, we can always normalise so that (-p)//2 < q <= p//2),
    then this sorting writes the entire collection of fibre parameters in a
    canonical form.
    """
    return tuple( sorted(fibreParams) )


def sortedFibreNegation(fibreParams):
    """
    Reverses the orientation of the given sequence of exceptional fibre
    parameters, and returns a sorted tuple of the reversed fibres.

    This routine assumes that the SFS is bounded, and always renormalises
    (2, q)-fibres as (2, 1)-fibres.
    """
    return sortedFibreParams(
            (2, 1) if p == 2 else (p, -q) for p, q in fibreParams )


def normalisedFibreParams(fibre):
    """
    Returns a normalised pair (p, q) for the given SFSFibre, such that
    (-p)//2 < q <= p//2.
    """
    p = fibre.alpha
    q = fibre.beta % p
    if q > p//2:
        q -= p
    return (p, q)


def fibrePreservingHomeomorphic( mySFS, yourSFS ):
    """
    Returns True if and only if the two given Seifert fibrations are
    fibre-preserving homeomorphic.

    The homeomorphism may be either orientation-preserving or
    orientation-reversing.

    Precondition:
    --> The two given Seifert fibrations are bounded and orientable.
    """
    if ( mySFS.baseClass() != yourSFS.baseClass() or
        mySFS.baseGenus() != yourSFS.baseGenus() or
        mySFS.punctures() != yourSFS.punctures() or
        mySFS.fibreCount() != yourSFS.fibreCount() ):
        return False

    # Same base surface and same numbers of fibres, so just need to compare
    # the fibre parameters. For this, we put everything in a canonical
    # normalised form, so that we can just do a direct equality check.
    myFibreParams = []
    yourFibreParams = []
    for i in range( mySFS.fibreCount() ):
        myFibreParams.append( normalisedFibreParams( mySFS.fibre(i) ) )
        yourFibreParams.append( normalisedFibreParams( yourSFS.fibre(i) ) )
    myFibreParams = sortedFibreParams(myFibreParams)
    return ( myFibreParams == sortedFibreParams(yourFibreParams) or
            myFibreParams == sortedFibreNegation(yourFibreParams) )
