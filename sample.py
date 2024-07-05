"""
Sample random knots from a table of knots.
"""
from sys import argv, stdout
from timeit import default_timer
from regina import *
from decomposeknot import decompose, DecompositionTracker


def extractFilenames(nameFile):
    """
    Extracts filenames from the given nameFile, which should be a text file
    with one filename per line.
    """
    sep = nameFile.rfind("/")
    if sep == -1:
        directory = ""
    else:
        directory = nameFile[:sep+1]
    output = []
    with open( nameFile, "r" ) as filenames:
        for line in filenames:
            output.append( directory + line.rstrip() )
    return output


def sample( size, *datasets ):
    """
    Returns a random sample of the given size from the knot tables in the
    given datasets.

    Each given dataset should be the name of a text file containing a list of
    filenames. Each listed filename should be a CSV file in the same
    directory as the corresponding dataset file, and each such csv file
    should include (at least) the following data:
    --> a column of knot names under the heading "name"; and
    --> a column of knot signatures under the heading "knot_sig".

    This routine returns a list of ( knot name, knot signature ) pairs,
    sampled uniformly at random, without replacement, from the given data.

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
    print( "-"*len(msg) )
    print()
    return sampleData


def decomposeFromSample( size, *datasets ):
    """
    Generates a random sample of the given size from the knot tables in the
    given datasets, and decomposes the knots in this sample.

    The knots are sampled using the sample() routine, so the format of the
    given datasets must adhere to the specifications stated in the
    documentation for sample().

    This routine prints the results of each decomposition to standard output.
    """
    strings = []
    for dataset in datasets:
        filenames = extractFilenames(dataset)
        for filename in filenames:
            strings.append(
                    filename.split( "/" )[-1].split( "." )[0] )
    title = ", ".join(strings)
    print()
    print( "+-" + "-"*len(title) + "-+" )
    print( "| {} |".format(title) )
    print( "+-" + "-"*len(title) + "-+" )
    print()
    stdout.flush()
    timedOutCases = []
    data = []
    totalTime = 0
    knots = sample( size, *datasets )
    for name, knotSig in knots:
        print(name)
        print( "-"*len(name) )
        knot = Link.fromKnotSig(knotSig)

        # Scale timeout time with the number of crossings.
        tracker = DecompositionTracker( True, 2*knot.size() )
        try:
            primes = decompose( knot, tracker )
        except TimeoutError as timeout:
            timedOutCases.append(name)
            print(timeout)
            print()
            continue
        if len(primes) == 0:
            print( "Unknot!" )
        elif len(primes) == 1:
            print( "Found 1 prime:" )
        else:
            print( "Found {} primes:".format( len(primes) ) )
        for i, loop in enumerate(primes):
            print( "    Drilled iso sig for prime #{}: {}".format(
                i, loop.drill().isoSig() ) )

        # Store data for post-processing.
        data.append( ( name, tracker.elapsed() ) )
        totalTime += tracker.elapsed()
        print()
    print( "="*32 )
    print( "Total knots: {}.".format( len(knots) ) )
    print( "Total time: {:.6f}.".format(totalTime) )
    if timedOutCases:
        print( "Cases that timed out ({} in total):".format(
            len(timedOutCases) ) )
        for name in timedOutCases:
            print( "    Name: {}.".format(name) )
    completedCount = len(knots) - len(timedOutCases)
    if completedCount:
        slowCoefficient = 2
        average = totalTime / completedCount
        print( "Cases slower than {} times the average:".format(
            slowCoefficient ) )
        noSlowCases = True
        for name, time in data:
            if time > slowCoefficient * average:
                noSlowCases = False
                print( "    Name: {}. Time: {:.6f}.".format( name, time ) )
        if noSlowCases:
            print( "    (None)" )
    print()
    return


def _sampleImpl( size, *datasets ):
    # Use reservoir sampling to limit memory usage to the size of the sample,
    # rather than the total size of all the given datasets.
    reservoir = []
    knotCount = 0
    for dataset in datasets:
        filenames = extractFilenames(dataset)
        for filename in filenames:
            with open( filename, "r" ) as table:
                headings = table.readline().rstrip().split( "," )
                nameCol = headings.index( "name" )
                sigCol = headings.index( "knot_sig" )
                while True:
                    row = table.readline()

                    # Does this row describe a knot?
                    if row == "":
                        # End of file.
                        break
                    entries = row.rstrip().split( "," )
                    knotSig = entries[sigCol]
                    if Link.fromKnotSig(knotSig).countComponents() != 1:
                        continue

                    # Reservoir sampling: the current knotSig should be
                    # included in the sample with probability
                    # (size // knotCount).
                    name = entries[nameCol]
                    knotCount += 1
                    if knotCount <= size:
                        reservoir.append( ( name, knotSig ) )
                    else:
                        i = RandomEngine.rand(knotCount)
                        if i < size:
                            reservoir[i] = ( name, knotSig )
    return ( reservoir, knotCount )


if __name__ == "__main__":
    decomposeFromSample( int( argv[1] ), *argv[2:] )
