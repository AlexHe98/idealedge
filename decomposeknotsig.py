"""
Decompose a knot given by a knot signature.
"""
from sys import argv
from regina import *
from decomposeknot import decompose


if __name__ == "__main__":
    # Run decompose() with the verbose option.
    decompose( Link.fromKnotSig( argv[1] ), True )
