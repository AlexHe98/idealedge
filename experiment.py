"""
Perform knot decomposition experiments in bulk.
"""
from sys import argv
from timeit import default_timer
from regina import *
from decomposeknot import decompose, DecompositionTracker


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


def decomposeFromTable( filename, skip=0, cap=None ):
    """
    Decomposes knots listed in the table in the given file, and prints the
    results to standard output.

    The given file should be a CSV file that includes (at least) the
    following data:
    --> a column of knot names under the heading "name"; and
    --> a column of knot signatures under the heading "knot_sig".

    It is possible to request that this routine skips the first few knots in
    the given table. By default, this routine does not skip any knots.

    This routine can be also run with an optional cap on the total number of
    knots that it processes. If the number of knots in the given table
    exceeds the cap, then this routine will terminate immediately after it
    has reached the cap.
    """
    title = filename.split( "/" )[-1].split( "." )[0]
    print()
    print( "+-" + "-"*len(title) + "-+" )
    print( "| {} |".format(title) )
    print( "+-" + "-"*len(title) + "-+" )
    print()
    knotCount = 0
    data = []
    totalTime = 0
    with open( filename, "r" ) as table:
        headings = table.readline().rstrip().split( "," )
        nameCol = headings.index( "name" )
        sigCol = headings.index( "knot_sig" )
        skipCount = 0
        while True:
            row = table.readline()
            if row == "":
                if knotCount == 0:
                    if skipCount > 0:
                        msg = "Skipped all knots."
                    else:
                        msg = "No knots found."
                    print( "-"*len(msg) )
                    print(msg)
                    print( "-"*len(msg) )
                    print()
                break
            entries = row.rstrip().split( "," )
            knot = Link.fromKnotSig( entries[sigCol] )
            if knot.countComponents() != 1:
                continue
            if skipCount < skip:
                skipCount += 1
                continue
            if skip > 0 and knotCount == 0:
                if skip == 1:
                    msg = "Skipped first knot."
                else:
                    msg = "Skipped first {} knots.".format(skip)
                print( "-"*len(msg) )
                print(msg)
                print( "-"*len(msg) )
                print()
            if cap is not None and knotCount >= cap:
                print( "----------------" )
                print( "Reached the cap." )
                print( "----------------" )
                print()
                break
            name = "Knot #{}: {}".format(
                    skip + knotCount, entries[nameCol] )
            print(name)
            print( "-"*len(name) )
            knotCount += 1
            tracker = DecompositionTracker()
            primes = decompose( knot, tracker )
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
    print( "Total knots: {}.".format(knotCount) )
    print( "Total time: {:.6f}.".format(totalTime) )
    if knotCount:
        slowCoefficient = 2
        average = totalTime / knotCount
        print( "Cases slower than {} times the average:".format(
            slowCoefficient ) )
        for name, time in data:
            if time > slowCoefficient * average:
                print( "    Name: {}. Time: {:.6f}.".format( name, time ) )
    print()
    return


if __name__ == "__main__":
    try:
        cap = int( argv[3] )
    except IndexError:
        cap = None
    try:
        skip = int( argv[2] )
    except IndexError:
        skip = 0
    decomposeFromTable( argv[1], skip, cap )
