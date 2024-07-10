"""
Perform knot decomposition experiments in bulk.
"""
from sys import argv, stdout
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


def readKnots(*datasets):
    """
    Reads knots from the knot tables in the given datasets.

    Each given dataset should be the name of a text file with one filename
    per line. Each listed filename should be a CSV file in the same directory
    as the corresponding dataset file. Each such CSV file should include (at
    least) the following data:
    --> a column of knot names under the heading "name"; and
    --> a column of knot signatures under the heading "knot_sig".
    For example, the knot tables available at
        https://regina-normal.github.io/data.html
    are given as CSV files that satisfy these requirements.

    For each knot in the given knot tables, this routine yields a pair
    (S, K), where:
    --> S is a string giving the name of the knot, as listed in the knot
        table; and
    --> K is the corresponding knot, encoded as a Regina Link object built
        directly from the knot signature listed in the knot table.
    """
    for dataset in datasets:
        for filename in extractFilenames(dataset):
            with open( filename, "r" ) as table:
                headings = table.readline().rstrip().split( "," )
                nameCol = headings.index( "name" )
                sigCol = headings.index( "knot_sig" )
                while True:
                    row = table.readline()
                    if row == "":
                        # End of file.
                        break

                    # Extract data from row.
                    entries = row.rstrip().split( "," )
                    knot = Link.fromKnotSig( entries[sigCol] )
                    if knot.countComponents() == 1:
                        yield ( entries[nameCol], knot )
    return


def runDecompositionExperiment(knotIterator):
    """
    Decomposes all knots described by the given iterator, and prints the
    results to standard output.

    This is a helper routine for running decomposition experiments in bulk.
    
    The given iterator should supply pairs of the form (S, K), where:
    --> S is a string giving a knot name; and
    --> K is a corresponding Regina Link object.
    """
    # Only want to keep the slow cases.
    slowCoefficient = 2
    slowKnots, slowTimes, timedOut, knotCount, totalTime = _experimentImpl(
            knotIterator, slowCoefficient )

    # Print summary.
    print( "="*32 )
    print( "Total knots: {}.".format(knotCount) )
    print( "Total time: {:.6f}.".format(totalTime) )
    if timedOut:
        print( "Cases that timed out ({} in total):".format(
            len(timedOut) ) )
        for name in timedOut:
            print(name)
    else:
        print( "All cases computed without timing out." )
    if slowKnots:
        print( "Cases slower than {} times the average:".format(
            slowCoefficient ) )
        for i, time, in enumerate(slowTimes):
            name = slowKnots[i][0]
            print( "    Time: {:.6f}. Name: {}.".format( time, name ) )
    else:
        print( "No cases slower than {} times the average.".format(
            slowCoefficient ) )
    print()

    # Rerun computation on slow knots (if any).
    if not slowKnots:
        return
    print( "="*32 )
    print( "Rerunning slow cases." )
    print()
    knots, times, timedOut, knotCount, totalTime = _experimentImpl(
            slowKnots, None )

    # Print summary.
    print( "="*32 )
    print( "Total knots: {}.".format(knotCount) )
    print( "Total time: {:.6f}.".format(totalTime) )
    if timedOut:
        print( "Cases that timed out ({} in total):".format(
            len(timedOut) ) )
        for name in timedOut:
            print(name)
    else:
        print( "All cases computed without timing out." )
    print()
    for i, oldTime in enumerate(slowTimes):
        newTime = times[i]
        name = knots[i][0]
        print( "Name: {}. Old time: {:.6f}. New time: {:.6f}.".format(
            name, oldTime, newTime ) )
    print()
    return


def _experimentImpl( knotIterator, slowCoefficient ):
    timedOutNames = []
    knots = []
    times = []
    for name, knot in knotIterator:
        knots.append( ( name, knot ) )
        print(name)
        print( "-"*len(name) )

        # Scale timeout time with the number of crossings.
        tracker = DecompositionTracker( True, knot.size() )
        try:
            primes = decompose( knot, tracker )
        except TimeoutError as timeout:
            timedOutNames.append(name)
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

        # Store elapsed time for post-processing.
        times.append( tracker.elapsed() )
        print()

    # Post-processing.
    knotCount = len(knots)
    totalTime = sum(times)
    if not slowCoefficient:
        return ( knots, times, timedOutNames, knotCount, totalTime )
    slowKnots = []
    slowTimes = []
    completedCount = knotCount - len(timedOutNames)
    if completedCount:
        average = totalTime / completedCount
        for i, time in enumerate(times):
            if time > slowCoefficient * average:
                slowKnots.append( knots[i] )
                slowTimes.append(time)
    return ( slowKnots, slowTimes, timedOutNames, knotCount, totalTime )


def decomposeKnots(*datasets):
    """
    Decomposes all knots from the knot tables in the given datasets, and
    prints the results to standard output.

    This routine uses the readKnots() routine to extract the knots from the
    given datasets, so the format of these datasets must adhere to the
    specifications stated in the documentation for readKnots().
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
    msg = "Decomposing knots from the following datasets:"
    print(msg)
    print( "="*len(msg) )
    printDatasetNames(*datasets)
    print()
    stdout.flush()
    runDecompositionExperiment( readKnots(*datasets) )
    return


def printDatasetNames(*datasets):
    """
    Prints the names of the given datasets to standard output.

    The format of the given datasets should adhere to the specifications
    stated in the documentation for the readKnots() routine.
    """
    for dataset in datasets:
        print( dataset.split( "/" )[-1].split( "." )[0] )
        for filename in extractFilenames(dataset):
            print( "    " + filename.split( "/" )[-1].split( "." )[0] )
    return


def decomposeAll( packet, verbose=True, insertChildren=False ):
    """
    Decomposes all knots appearing as descendants of the given packet.

    With the verbose option (which is switched on by default), this routine
    prints regular reports on the progress of the experiments. With the
    insertChildren option (which is switched *off* by default), this routine
    inserts the results of each decomposition as a child of the corresponding
    PacketOfLink object.
    """
    for knot in packet.descendants():
        if ( not isinstance( knot, PacketOfLink ) or
                knot.countComponents() != 1 ):
            continue
        decompose( knot, verbose, insertChildren )
        if verbose:
            print()
    return


def clearExperiments(packet):
    """
    Deletes all descendants of the given packet whose labels begin with the
    string "Primes".
    """
    doomed = []
    for child in packet.children():
        clearExperiments(child)
        if child.label()[:6] == "Primes":
            doomed.append(child)
    for d in doomed:
        d.makeOrphan()
    return


if __name__ == "__main__":
    decomposeKnots( *argv[1:] )
