"""
Run recogniseSFS() on triangulations built by gluing a SFS to a hyperbolic
3-manifold.

Usage: Must supply a path to which to write the experimental data. To avoid
    overwriting data, this script will check that the path is to a file which
    does not yet exist.
"""
from pathlib import Path
import sys
from timeit import default_timer
from regina import *
from recsfs import recogniseSFS, SFSRecognitionTracker


if __name__ == "__main__":
    outputPath = sys.argv[1]
    useHeuristics = True
    with ( open( outputPath, 'x' ) as outputFile,
          open( "data/nonsfs/sample.txt", 'r' ) as sigLines ):
        for line in sigLines:
            neoSig = line.rstrip()
            tri = Triangulation3.fromSig(neoSig)
            size = tri.size()
            tracker = SFSRecognitionTracker()
            start = default_timer()
            ans = recogniseSFS( tri, useHeuristics, tracker )
            time = default_timer() - start
            msg = "{} {} {:.6f}".format(
                    size,
                    int( tracker.usedCombinatorialRecognition() ),
                    time )
            outputFile.write( msg + "\n" )
            outputFile.flush()
            print(msg)
            sys.stdout.flush()

            # Check that we actually got the correct answer.
            assert ans is None
