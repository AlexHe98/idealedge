"""
Generate random composite knots by sampling from knot tables.
"""
from sys import argv, stdout
from regina import *
from sample import sample
from decomposeknot import decompose, DecompositionTracker


def generateComposites( numKnots, numSummands, *filenames ):
    """
    Generates the given number of knots, each given by composing the given
    number of summands, by randomly sampling summands from the knot table(s)
    in the given file(s).

    Each given file should be a CSV file that includes (at least) the
    following data:
    --> a column of knot names under the heading "name"; and
    --> a column of knot signatures under the heading "knot_sig".
    """
    knotNamesAndSigs = sample( numKnots * numSummands, *filenames )
    numSampled = len(knotNamesAndSigs)

    # Before generating composites, use Fisher-Yates to shuffle the sample.
    for i in range( numSampled - 1, 0, -1 ):
        swapi = RandomEngine.rand(i+1)
        knotNamesAndSigs[i], knotNamesAndSigs[swapi] = (
                knotNamesAndSigs[swapi], knotNamesAndSigs[i] )

    # Use the shuffled summands to generate a list of composite knots.
    numKnots = min( numKnots, numSampled // numSummands )
    composites = []
    for i in range(numKnots):
        compositeName, knotSig = knotNamesAndSigs[ i*numSummands ]
        compositeKnot = Link.fromKnotSig(knotSig)
        for j in range( i*numSummands + 1, (i+1)*numSummands ):
            name, knotSig = knotNamesAndSigs[j]
            compositeName += " # {}".format(name)
            compositeKnot.composeWith( Link.fromKnotSig(knotSig) )
        composites.append( ( compositeName, compositeKnot ) )
    return composites


def decomposeComposites( numKnots, numSummands, *filenames ):
    """
    Generates the given number of knots, each given by composing the given
    number of summands, and then decomposes all the generated knots.

    The summands are sampled from the knot table(s) in the given files(s).
    Each such file should be a CSV file that includes (at least) the
    following data:
    --> a column of knot names under the heading "name"; and
    --> a column of knot signatures under the heading "knot_sig".

    The routine prints the results of each decomposition to standard output.
    """
    title = "Knots given by composing {} summands".format(numSummands)
    print()
    print( "+-" + "-"*len(title) + "-+" )
    print( "| {} |".format(title) )
    print( "+-" + "-"*len(title) + "-+" )
    print()
    print( "Summands sampled from the following data sets:" )
    for filename in filenames:
        print( "    {}".format(
            filename.split( "/" )[-1].split( "." )[0] ) )
    print()
    stdout.flush()
    timedOutCases = []
    data = []
    totalTime = 0
    knots = generateComposites( numKnots, numSummands, *filenames )
    for name, knot in knots:
        print(name)
        print( "-"*len(name) )

        # Scale timeout time with the number of crossings.
        tracker = DecompositionTracker( True, knot.size() )
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
    print( "Total knots: {}.".format(numKnots) )
    print( "Total time: {:.6f}.".format(totalTime) )
    if timedOutCases:
        print( "Cases that timed out ({} in total):".format(
            len(timedOutCases) ) )
        for name in timedOutCases:
            print( "    Name: {}.".format(name) )
    completedCount = numKnots - len(timedOutCases)
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


if __name__ == "__main__":
    numKnots = int( argv[1] )
    numSummands = int( argv[2] )
    decomposeComposites( numKnots, numSummands, *argv[3:] )
