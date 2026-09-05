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


def hardSFSPath( isBaseOrbl=None, genus=None, numBdries=None ):
    """
    Returns the path to a file containing non-standard triangulations of
    bounded orientable Seifert fibred spaces with the given properties.

    If such a file exists, then it should have been generated
    programmatically, and each line of the file should follow the format
    described in sfsLine().

    If no properties are provided, then this routine simply returns the path
    to the directory containing all such files.
    """
    if isBaseOrbl is None:
        return "sfs-samples/"
    elif isBaseOrbl:
        base = "baseOr"
    else:
        base = "baseNor"
    if (genus is None) or (numBdries is None):
        raise ValueError( "hardSFSPath() requires that either all three " +
                         "properties are supplied, or none at all" )
    return "sfs-samples/sfs-{}-g{}-m{}.txt".format( base, genus, numBdries )


def parseSFSFileName(fileName):
    """
    Parses a file name of the format 'sfs-<base>-g<genus>-m<numBdries>.txt'.

    This routine returns a tuple (isBaseOrbl, genus, numBdries), where:
    --> isBaseOrbl is True if <base> is "baseOr", and False if <base> is
        "baseNor";
    --> genus is a non-negative (moreover, positive if isBaseOrbl is False)
        integer; and
    --> numBdries is a positive integer.
    """
    if fileName[-4:] != ".txt":
        raise ValueError(
                "parseSFSFileName() requires '.txt' file extension" )
    split = fileName[:-4].split("-")
    msg = "Incorrect format for parseSFSFileName()"
    if ( len(split) != 4 or split[0] != "sfs" or split[1][:4] != "base" or
        split[2][0] != "g" or split[3][0] != "m" ):
        raise ValueError(msg)
    try:
        # Because we split at "-" characters, any integers we get should be
        # non-negative.
        genus = int( split[2][1:] )
        numBdries = int( split[3][1:] )
    except Exception:
        raise ValueError(msg)
    if split[1][4:] == "Or":
        isBaseOrbl = True
    elif split[1][4:] == "Nor":
        isBaseOrbl = False
        if genus == 0:
            raise ValueError(
                    "Non-orientable base surface cannot have genus 0" )
    else:
        raise ValueError(msg)
    if numBdries == 0:
        raise ValueError( "Number of boundaries should be positive" )
    return (isBaseOrbl, genus, numBdries)
