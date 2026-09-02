"""
Scratch work for experimental testing of recogniseSFS().
"""
import sys
from math import gcd as pythonGCD
from itertools import combinations_with_replacement as combinsWithRep


if __name__ == "__main__":
#    baseSignedGenus = int( sys.argv[1] )
#    numBdries = int( sys.argv[2] )
#    numFibres = int( sys.argv[3] )
    maxMultiplicity = 5 #TODO
    numFibres = int( sys.argv[1] )  #TODO
    fibres = []
    for p in range( 2, 1 + maxMultiplicity ):
        for q in range( 1, p ):
            if pythonGCD(p, q) != 1:
                continue
            fibres.append( (p, q) )
    count = 0
    for res in combinsWithRep( fibres, numFibres ):
        count += 1
        print( count, res )
