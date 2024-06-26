"""
Sample random knots from a table of knots.
"""
from sys import argv
from timeit import default_timer
from regina import *
from decomposeknot import decompose, DecompositionTracker


def sample( size=1, *filenames ):
    """
    Returns a random sample of the given size from the knot table(s) in the
    given file(s).

    Each given file should be a CSV file that includes (at least) the
    following data:
    --> a column of knot names under the heading "name"; and
    --> a column of knot signatures under the heading "knot_sig".
    This routine returns a list of ( knot name, knot signature ) pairs,
    sampled uniformly at random, without replacement, from the given data.
    """
    start = default_timer()
    RandomEngine.reseedWithHardware()
    reservoir = []
    knotCount = 0
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

                # Reservoir sampling: ensures that the current knotSig is
                # included in the sample with probability (size//knotCount).
                name = entries[nameCol]
                knotCount += 1
                if knotCount <= size:
                    reservoir.append( ( name, knotSig ) )
                else:
                    i = RandomEngine.rand(knotCount)
                    if i < size:
                        reservoir[i] = ( name, knotSig )

    # Finished sampling.
    msg = "Sampled {} out of {} knots. Time: {:.6f}.".format(
        min( size, knotCount ), knotCount, default_timer() - start )
    print(msg)
    print( "-"*len(msg) )
    print()
    return reservoir


def decomposeFromSample( size, *filenames ):
    """
    Generates a random sample of the given size from the knot table(s) in the
    given file(s), and decomposes the knots in this sample.

    Each given file should be a CSV file that includes (at least) the
    following data:
    --> a column of knot names under the heading "name"; and
    --> a column of knot signatures under the heading "knot_sig".

    This routine prints the results of each decomposition to standard output.
    """
    title = ""
    for filename in filenames:
        title += ", " + filename.split( "/" )[-1].split( "." )[0]
    title = title[2:]
    print()
    print( "+-" + "-"*len(title) + "-+" )
    print( "| {} |".format(title) )
    print( "+-" + "-"*len(title) + "-+" )
    print()
    data = []
    totalTime = 0
    knots = sample( size, *filenames )
    for name, knotSig in knots:
        print(name)
        print( "-"*len(name) )
        tracker = DecompositionTracker(True)
        primes = decompose( Link.fromKnotSig(knotSig), tracker )
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
    if knots:
        slowCoefficient = 2
        average = totalTime / len(knots)
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


if __name__ == "__main__":
    decomposeFromSample( int( argv[1] ), *argv[2:] )
