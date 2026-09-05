"""
Scratch work for SFS experiments.
"""
from timeit import default_timer
from regina import *
from recsfs import recogniseSFS
from experiments.sfs.io import hardSFSPath, readSFSLine


if __name__ == "__main__":
    #TODO
    isBaseOrbl = False
    baseGenus = 4
    numBdries = 4
    filePath = hardSFSPath( isBaseOrbl, baseGenus, numBdries )
    with open( filePath, 'r' ) as sfsFile:
        for line in sfsFile:
            neoSig, fibres = readSFSLine(line)
            tri = Triangulation3.fromSig(neoSig)
            size = tri.size()
            #TODO Compare computed vs expected results.
            start = default_timer()
            ans = recogniseSFS(tri)
            time = default_timer() - start
            print( "Size: {}. Time: {:.6f}. {}".format( size, time, ans ) )
