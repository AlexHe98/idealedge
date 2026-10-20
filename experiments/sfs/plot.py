"""
Plot results from SFS experiments.
"""
from matplotlib import pyplot
from pandas import DataFrame
from seaborn import scatterplot
import sys


if __name__ == "__main__":
    dataPath = sys.argv[1]
    try:
        yscale = sys.argv[2]
    except IndexError:
        yscale = "linear"
    try:
        colourCombRec = bool( int( sys.argv[3] ) )
    except IndexError:
        colourCombRec = False
    sizes = []
    times = []
    usedCombRec = []
    with open( dataPath, 'r' ) as dataFile:
        for line in dataFile:
            size, combRecIndicator, time = line.rstrip().split(" ")
            sizes.append( int(size) )
            usedCombRec.append( bool( int(combRecIndicator) ) )
            times.append( float(time) )
    sizeName = "Number of tetrahedra"
    timeName = "Time (seconds)"
    combRecName = "Early termination using combinatorial recognition"
    data = DataFrame( {
        sizeName: sizes,
        timeName: times,
        "combRecRaw": usedCombRec } )
    combRecCounts = data["combRecRaw"].value_counts()
    data[combRecName] = data["combRecRaw"].map(
            lambda combRecBool:
            f"{combRecBool} (n={combRecCounts[combRecBool]})" )
    _, ax = pyplot.subplots( figsize=(8, 4.5) )
    if colourCombRec:
        scatterplot( data=data, x=sizeName, y=timeName, ax=ax,
                    hue=combRecName )
    else:
        scatterplot( data=data, x=sizeName, y=timeName, ax=ax )
    pyplot.title("Timings for bounded orientable SFS recognition")
    pyplot.yscale(yscale)
    pyplot.show()
