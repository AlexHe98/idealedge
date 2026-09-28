"""
Scratch work for SFS experiments.
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
    scatterplot( data=data, x=sizeName, y=timeName,
                hue=altEnumName, palette="flare_r" )
    pyplot.title("Timings for bounded orientable SFS recognition")
    pyplot.yscale(yscale)
    pyplot.show()
