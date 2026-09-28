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
        useAltEnum = bool( int( sys.argv[3] ) )
    except IndexError:
        useAltEnum = False
    sizes = []
    times = []
    altEnumCounts = []
    with open( dataPath, 'r' ) as dataFile:
        for line in dataFile:
            size, altEnums, time = line.rstrip().split(" ")
            sizes.append( int(size) )
            altEnumCounts.append( int(altEnums) )
            times.append( float(time) )
    sizeName = "Number of tetrahedra"
    timeName = "Time (seconds)"
    altEnumName = "Number of alternate enumerations used"
    data = DataFrame( {
        sizeName: sizes,
        timeName: times,
        altEnumName: altEnumCounts } )
    if useAltEnum:
        scatterplot( data=data, x=sizeName, y=timeName,
                    hue=altEnumName, palette="flare_r" )
    else:
        scatterplot( data=data, x=sizeName, y=timeName )
    pyplot.title("Timings for bounded orientable SFS recognition")
    pyplot.yscale(yscale)
    pyplot.show()
