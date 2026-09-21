"""
Run recogniseSFS() on hard triangulations of bounded orientable SFS.

Here, "hard" means that combinatorial recognition heuristics fail to identify
the Seifert fibred structure, which makes such hard triangulations good test
cases for the practical performance of the bounded orientable SFS recognition
algorithm implemented by the recogniseSFS() routine.

Usage: Must supply a path to which to write the experimental data. To avoid
    overwriting data, this script will check that the path is to a file which
    does not yet exist.
"""
from pathlib import Path
import sys
from timeit import default_timer
from regina import *
from aux.sfs import fibrePreservingHomeomorphic
from recsfs import recogniseSFS, SFSRecognitionTracker
from experiments.sfs.io import hardSFSPath, readSFSLine, parseSFSFileName


if __name__ == "__main__":
    outputPath = sys.argv[1]
    useHeuristics = True
    with open( outputPath, 'x' ) as outputFile:
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
                        msg = "{} {} {:.6f}".format(
                                size, tracker.alternateEnumerationsCount(),
                                time )
                        outputFile.write( msg + "\n" )
                        outputFile.flush()
                        print(msg)
                        sys.stdout.flush()

                        # Check that we actually got the correct answer.
                        expected = SFSpace( baseClass, genus, numBdries )
                        for p, q in fibres:
                            expected.insertFibre( p, q )
                        assert fibrePreservingHomeomorphic( expected, ans ),\
                                expected
