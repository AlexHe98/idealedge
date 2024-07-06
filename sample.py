"""
Sample random knots from a table of knots, and perform knot decomposition
experiments on such samples.
"""
from sys import argv, stdout
from timeit import default_timer
from regina import *
from experiment import runDecompositionExperiment
from experiment import readKnots, printDatasetNames


def sample( size, *datasets ):
    """
    Returns a random sample of the given size from the knot tables in the
    given datasets.

    This routine uses the readKnots() routine to extract the knots from the
    given datasets, so the format of these datasets must adhere to the
    specifications stated in the documentation for readKnots().

    This routine returns a list of pairs of the form (S, K), where:
    --> S is a string giving a knot name; and
    --> K is a corresponding Regina Link object.

    If the requested sample size is larger than the total number of knots in
    the given data, then this routine will simply repeatedly sample the
    entire data set until it reaches the required sample size.
    """
    if size < 1:
        raise ValueError( "Sample size must be positive." )
    start = default_timer()
    RandomEngine.reseedWithHardware()
    sampleData, knotCount = _sampleImpl( size, *datasets )
    if knotCount == 0:
        # Just in case the given datasets do not contain any knots.
        return []

    # Keep sampling until we have enough knots.
    remaining = size - knotCount
    while remaining > 0:
        sampleData.extend( _sampleImpl( remaining, *datasets )[0] )
        remaining -= knotCount
    msg = "Sampled {} knots from set of {}. Time: {:.6f}.".format(
            size, knotCount, default_timer() - start )
    print(msg)
    print( "="*len(msg) )
    print()
    return sampleData


def decomposeFromSample( size, *datasets ):
    """
    Generates a random sample of the given size from the knot tables in the
    given datasets, and decomposes the knots in this sample.

    The knots are sampled using the sample() routine. This means that the
    format of the given datasets must adhere to the specifications stated in
    the documentation for the readKnots() routine.

    This routine prints the results of each decomposition to standard output.
    """
    strings = []
    for dataset in datasets:
        strings.append(
                dataset.split( "/" )[-1].split( "." )[0] )
    title = ", ".join(strings)
    print()
    print( "+-" + "-"*len(title) + "-+" )
    print( "| {} |".format(title) )
    print( "+-" + "-"*len(title) + "-+" )
    print()
    msg = "Decomposing a sample of knots from the following datasets:"
    print(msg)
    print( "="*len(msg) )
    printDatasetNames(*datasets)
    print()
    stdout.flush()
    runDecompositionExperiment( sample( size, *datasets ) )
    return


def _sampleImpl( size, *datasets ):
    # Use reservoir sampling to limit memory usage to the size of the sample,
    # rather than the total size of all the given datasets.
    reservoir = []
    knotCount = 0
    for name, knot in readKnots(*datasets):
        knotCount += 1
        if knotCount <= size:
            reservoir.append( ( name, knot ) )
        else:
            i = RandomEngine.rand(knotCount)
            if i < size:
                reservoir[i] = ( name, knot )
    return ( reservoir, knotCount )


if __name__ == "__main__":
    decomposeFromSample( int( argv[1] ), *argv[2:] )
