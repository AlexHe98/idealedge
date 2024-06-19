"""
Perform knot decomposition experiments in bulk.
"""
from sys import argv
from timeit import default_timer
from regina import *
from decomposeknot import decompose


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


def decomposeAllInTable( filename, verbose=True ):
    """
    Decomposes all knots listed in the table in the given file.

    The given file should be a CSV file that includes (at least) the
    following data:
    --> a column of knot names under the heading "name"; and
    --> a column of knot signatures under the heading "knot_sig".

    With the verbose option (which is switched on by default), this routine
    prints regular reports on the progress of the experiments.
    """
    if verbose:
        name = filename.split( "/" )[-1].split( "." )[0]
        print(name)
        print( "="*len(name) )
        print()
        knotCount = 0
        start = default_timer()
    with open( filename, "r" ) as table:
        headings = table.readline().rstrip().split( "," )
        nameCol = headings.index( "name" )
        sigCol = headings.index( "knot_sig" )
        while True:
            row = table.readline()
            if row == "":
                break
            entries = row.rstrip().split( "," )
            knot = Link.fromKnotSig( entries[sigCol] )
            if knot.countComponents() != 1:
                continue
            if verbose:
                name = "Knot #{}: {}".format(
                        knotCount, entries[nameCol] )
                print(name)
                print( "-"*len(name) )
                knotCount += 1
            primes = decompose( knot, verbose )
            if verbose:
                if len(primes) == 0:
                    print( "Unknot!" )
                elif len(primes) == 1:
                    print( "Found 1 prime:" )
                else:
                    print( "Found {} primes:".format( len(primes) ) )
                for i, loop in enumerate(primes):
                    print( "    Drilled iso sig for prime #{}: {}".format(
                        i, loop.drill().isoSig() ) )
                print()
    print( "="*32 )
    print( "Total knots: {}.".format(knotCount) )
    print( "Total time: {:.6f}.".format( default_timer() - start ) )
    return


if __name__ == "__main__":
    decomposeAllInTable( argv[1] )
