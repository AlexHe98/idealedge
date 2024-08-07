"""
Generate random composite knots by sampling summands from knot tables, and
perform knot decomposition experiments on the knots generated in this way.
"""
from sys import argv, stdout
from regina import *
from experiment import runDecompositionExperiment
from experiment import printDatasetNames
from sample import sample


def generateComposites( numKnots, numSummands, *datasets ):
    """
    Generates the given number of knots, each given by composing the given
    number of summands, by sampling summands from the given datasets.

    Each summand is generated by selecting one of the datasets uniformly at
    random, and then choosing one of the knots in the selected dataset
    uniformly at random.

    This routine uses the sample() routine to sample from each individual
    dataset. This means that the format of the given datasets must adhere to
    the specifications stated in the documentation for the readKnots()
    routine.

    For each generated composite knot, this routine yields a pair of the form
    (S, K), where:
    --> S is a string giving a knot name; and
    --> K is a corresponding Regina Link object.
    """
    numSets = len(datasets)
    if numSets == 1:
        sampleSizes = [ numKnots * numSummands ]
        summandSources = [ [0]*numSummands for _ in range(numKnots) ]
    else:
        RandomEngine.reseedWithHardware()

        # Which dataset should we get each summand from?
        sampleSizes = [0]*numSets
        summandSources = []
        for _ in range(numKnots):
            temp = []
            for _ in range(numSummands):
                # Randomly pick one of the datasets from which to generate
                # one of the summands.
                setNumber = RandomEngine.rand(numSets)
                sampleSizes[setNumber] += 1
                temp.append(setNumber)
            summandSources.append(temp)

    # Randomly sample summands from each dataset.
    summandNamesAndKnots = []
    for setNumber in range(numSets):
        summandSample = sample(
                sampleSizes[setNumber], datasets[setNumber] )
        summandNamesAndKnots.append(summandSample)

        # Use Fisher-Yates to shuffle the sample (in-place).
        for i in range( sampleSizes[setNumber] - 1, 0, -1 ):
            swapi = RandomEngine.rand(i+1)
            summandSample[i], summandSample[swapi] = (
                    summandSample[swapi], summandSample[i] )

    # Use sampled summands to generate a list of composite knots.
    composites = []
    for sources in summandSources:
        summandNames = []
        summandKnots = []
        for setNumber in sources:
            name, knot = summandNamesAndKnots[setNumber].pop()
            summandNames.append(name)
            summandKnots.append(knot)

        # Build the composite knot.
        compositeName = summandNames[0]
        compositeKnot = summandKnots[0]
        for i in range( 1, numSummands ):
            compositeName += " # {}".format( summandNames[i] )
            compositeKnot.composeWith( summandKnots[i] )
        composites.append( ( compositeName, compositeKnot ) )
    return composites


def decomposeComposites( numKnots, numSummands, *datasets ):
    """
    Generates the given number of knots, each given by composing the given
    number of summands, and then decomposes all the generated knots.

    The composite knots are generated using the generateComposites() routine.
    This means that the format of the given datasets must adhere to the
    specifications stated in the documentation for the readKnots() routine.

    This routine prints the results of each decomposition to standard output.
    """
    title = "Knots given by composing {} summands".format(numSummands)
    print()
    print( "+-" + "-"*len(title) + "-+" )
    print( "| {} |".format(title) )
    print( "+-" + "-"*len(title) + "-+" )
    print()
    msg = "Sampling summands from the following data sets:"
    print(msg)
    print( "="*len(msg) )
    printDatasetNames(*datasets)
    print()
    stdout.flush()
    runDecompositionExperiment( generateComposites(
        numKnots, numSummands, *datasets ) )
    return


if __name__ == "__main__":
    numKnots = int( argv[1] )
    numSummands = int( argv[2] )
    decomposeComposites( numKnots, numSummands, *argv[3:] )
