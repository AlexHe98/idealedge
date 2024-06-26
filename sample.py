"""
Sample random knots from a table of knots.
"""
from sys import argv
from timeit import default_timer
from regina import *
from decomposeknot import decompose, DecompositionTracker


def sample( filename, size=1 ):
    """
    Returns a random sample of the given size from the knot table in the
    given file.

    The given file should be a CSV file that includes (at least) the
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
    with open( filename, "r" ) as table:
        headings = table.readline().rstrip().split( "," )
        nameCol = headings.index( "name" )
        sigCol = headings.index( "knot_sig" )
        while True:
            row = table.readline()

            # Does this row describe a knot?
            if row == "":
                # End of file.
                msg = "Sampled {} out of {} knots. Time: {:.6f}.".format(
                    min( size, knotCount ), knotCount,
                    default_timer() - start )
                print(msg)
                print( "-"*len(msg) )
                print()
                return reservoir
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


def decomposeFromSample( filename, size ):
    """
    Generates a random sample of the given size from the knot table in the
    given file, and decomposes the knots in this sample.

    The given file should be a CSV file that includes (at least) the
    following data:
    --> a column of knot names under the heading "name"; and
    --> a column of knot signatures under the heading "knot_sig".

    This routine prints the results of each decomposition to standard output.
    """
    title = filename.split( "/" )[-1].split( "." )[0]
    print()
    print( "+-" + "-"*len(title) + "-+" )
    print( "| {} |".format(title) )
    print( "+-" + "-"*len(title) + "-+" )
    print()
    data = []
    totalTime = 0
    knots = sample( filename, size )
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
        for name, time in data:
            if time > slowCoefficient * average:
                print( "    Name: {}. Time: {:.6f}.".format( name, time ) )
    print()
    return


if __name__ == "__main__":
    try:
        size = int( argv[2] )
    except IndexError:
        size = 1
    decomposeFromSample( argv[1], size )
