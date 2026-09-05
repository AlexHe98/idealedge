"""
Scratch work for SFS experiments.
"""
from pathlib import Path
from timeit import default_timer
from regina import *
from recsfs import recogniseSFS
from experiments.sfs.io import hardSFSPath, readSFSLine, parseSFSFileName


if __name__ == "__main__":
    for file in Path( hardSFSPath() ).iterdir():
        if file.is_file():
            isBaseOrbl, genus, numBdries = parseSFSFileName(file.name)
            with open( hardSFSPath(
                isBaseOrbl, genus, numBdries ) ) as sfsFile:
                for line in sfsFile:
                    neoSig, fibres = readSFSLine(line)
                    tri = Triangulation3.fromSig(neoSig)
                    size = tri.size()
                    #TODO Compare computed vs expected results.
                    start = default_timer()
                    ans = recogniseSFS(tri)
                    time = default_timer() - start
                    print( "Size: {}. Time: {:.6f}. {}".format(
                        size, time, ans ) )
            #TODO Temporarily break out while testing.
            break
