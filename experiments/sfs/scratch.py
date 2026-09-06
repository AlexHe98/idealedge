"""
Scratch work for SFS experiments.
"""
from pathlib import Path
from timeit import default_timer
from regina import *
from aux.sfs import fibrePreservingHomeomorphic
from recsfs import recogniseSFS, SFSRecognitionTracker
from experiments.sfs.io import hardSFSPath, readSFSLine, parseSFSFileName


if __name__ == "__main__":
    useHeuristics = True
    for file in Path( hardSFSPath() ).iterdir():
        if file.is_file():
            isBaseOrbl, genus, numBdries = parseSFSFileName(file.name)
            if isBaseOrbl:
                baseClass = SFSpace.Class.bo1
            else:
                baseClass = SFSpace.Class.bn2
            with open( hardSFSPath(
                isBaseOrbl, genus, numBdries ) ) as sfsFile:
                for line in sfsFile:
                    neoSig, fibres = readSFSLine(line)
                    tri = Triangulation3.fromSig(neoSig)
                    size = tri.size()
                    tracker = SFSRecognitionTracker()
                    start = default_timer()
                    ans = recogniseSFS( tri, useHeuristics, tracker )
                    time = default_timer() - start
                    print( "Size: {}. Time: {:.6f}. #altEnums: {}. {}".format(
                        size, time,
                        tracker.alternateEnumerationsCount(), ans ) )

                    # Check that we actually got the correct answer.
                    expected = SFSpace( baseClass, genus, numBdries )
                    for p, q in fibres:
                        expected.insertFibre( p, q )
                    assert fibrePreservingHomeomorphic( expected, ans ),\
                            expected
            #TODO Temporarily break out while testing.
            break
