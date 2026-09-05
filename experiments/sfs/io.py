"""
File I/O for bounded orientable SFS experiments.
"""


def sfsLine( neoSig, fibres ):
    """
    Returns a string containing the given data.

    The given neoSig should be a Regina 2nd-generation iso sig. Each element
    of fibres (which is allowed to be empty) should be a pair (p, q) of
    integers specifying valid parameters for an exceptional fibre in a
    Seifert fibred space (i.e., p >= 2 and gcd(p, q) == 1). The returned
    string will end with a single newline, and will be a space-separated list
    containing the neoSig followed by each fibre in the format "(p,q)" (with
    *no space* after the comma).
    """
    line = neoSig
    for p, q in fibres:
        line += " ({},{})".format( p, q )
    return line + "\n"


def readSFSLine(line):
    """
    Reads a 2nd-generation iso sig and a tuple of fibre parameters from the
    given line.

    The line should follow the format described in sfsLine().
    """
    data = line.rstrip().split(" ")
    fibres = []
    neoSig = data[0]
    for fibreString in data[1:]:
        fibres.append( tuple(
            int(k) for k in fibreString[1:-1].split(",") ) )
    return neoSig, tuple(fibres)


def hardSFSPath( isBaseOrbl, genus, numBdries ):
    """
    Returns the path to a file containing non-standard triangulations of
    bounded orientable Seifert fibred spaces with the given properties.

    If such a file exists, then it should have been generated
    programmatically, and each line of the file should follow the format
    described in sfsLine().
    """
    if isBaseOrbl:
        base = "baseOr"
    else:
        base = "baseNor"
    return "sfs-samples/sfs-{}-g{}-m{}.txt".format( base, genus, numBdries )
