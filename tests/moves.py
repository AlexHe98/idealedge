"""
Test suite for 2-3, 3-2 and 2-0 moves.
"""
from regina import *
from edgelabel import EdgeLabelling
from moves import twoThree, threeTwo, twoZero, fourFour, twoOne
from sys import argv, stdout
from tests.aux import parseTestNames, allTestsPassedMessage
from tests.aux import doTest, runNamedTestSuite

#NOTE These tests use the Isomorphism3 bracket operator, introduced in
#       Regina 7.1, to apply an isomorphism to a Triangulation3 object. In
#       older versions of Regina, equivalent functionality was provided by the
#       Isomorphism3.apply() routine.


def _generateIsomorphisms( triSize, maxIsos ):
    iso = Isomorphism3.identity(triSize)
    iso.inc()
    count = 0
    while ( count < maxIsos ) and ( not iso.isIdentity() ):
        count += 1
        yield iso
        for _ in range(count*25):
            iso.inc()
    return

def _multiplicities(edge):
    """
    Returns a list of the multiplicities of the given edge, sorted from
    largest to smallest.

    The multiplicity of an edge e with respect to a tetrahedron T is the
    number of times (max 6) in which e is embedded as an edge of T. The
    returned list will contain all nonzero multiplicities of the given edge.

    Note that the sum of the multiplicities is precisely the degree.
    """
    tri = edge.triangulation()
    ans = []
    for tet in tri.tetrahedra():
        mult = 0
        for e in range(6):
            if tet.edge(e) == edge:
                mult += 1
        if mult:
            ans.append(mult)
    ans.sort( reverse=True )
    return ans

def _relativeEdgeOrientations( tri, iso ):
    """
    Returns the relative edge orientations resulting from applying the given
    isomorphism to the given triangulation.
    """
    newTri = iso(tri)
    ans = []
    for e in tri.edges():
        tetIndex = e.front().tetrahedron().index()
        oldVerts = e.front().vertices()
        tetImage = iso.simpImage(tetIndex)
        newVerts = iso.facetPerm(tetIndex) * oldVerts

        # Compare newVerts with the vertex permutation given by edgeMapping().
        compareVerts = newTri.tetrahedron(tetImage).edgeMapping(
                Edge3.faceNumber(newVerts) )
        if compareVerts[0] == newVerts[0]:
            # Same orientation.
            ans.append(1)
        else:
            # Opposite orientation.
            ans.append(-1)
    return ans

def _verifyRelab( before, relab, inter, inlab, after ):
    """
    Checks that relab and inlab appropriately track edge relabellings
    resulting respectively from a move taking before to inter, and then a
    supposed inverse move taking inter to after.
    """
    # The idea is to search for the isomorphism that suitably relates the
    # before and after triangulations.
    moveRelOr = []
    for i in range( before.countEdges() ):
        # Relative orientation after initial move.
        tetIndex = relab[i].tetrahedron().index()
        verts = relab[i].vertices()
        compareVerts = inter.tetrahedron(tetIndex).edgeMapping(
                Edge3.faceNumber(verts) )
        if verts[0] == compareVerts[0]:
            moveRelOr.append(1)
        else:
            moveRelOr.append(-1)

        # Relative orientation after inverse move.
        ii = relab.underlyingEdgeIndex(i)
        tetIndex = inlab[ii].tetrahedron().index()
        verts = inlab[ii].vertices()
        compareVerts = after.tetrahedron(tetIndex).edgeMapping(
                Edge3.faceNumber(verts) )
        if verts[0] != compareVerts[0]:
            moveRelOr[-1] *= -1
    for iso in before.findAllIsomorphisms(after):
        if moveRelOr == _relativeEdgeOrientations( before, iso ):
            # Found the desired isomorphism.
            return True
    return False


def testTwoThree():
    total = 0
    for testSig in [ "cMcabbgqs",               # 2 2-3 moves
                    "gLLPQcdefeffpvauppb" ]:    # 12 2-3 moves
        total += _test23all(testSig)
    return total

def _test23all( testSig, maxIsos=16 ):
    print( "2-3 and 3-2 moves on \"{}\"".format(testSig) )
    stdout.flush()
    t = Triangulation3.fromIsoSig(testSig)
    t.orient()

    # Test 2-3 moves on all eligible triangles.
    count = 0
    for f in t.triangles():
        # Make sure that we get the correct isomorphism type, and that we
        # preserve orientation.
        #
        #NOTE Triangulation3.withPachner(f) was introduced in Regina 7.4.
        #       Older versions of Regina did not provide a routine to
        #       construct a new triangulation via a 2-3 move; instead, this
        #       behaviour was achieved by first constructing a copy of the
        #       triangulation, and then performing the appropriate 2-3 move on
        #       the copy.
        pach = t.withPachner(f)
        if pach is None:
            # 2-3 move not eligible on f.
            continue
        count += 1
        try:
            _test23single( f, pach )
        except AssertionError as ae:
            print(t)
            raise ae

        # To test as many cases of the implementation as possible, test the
        # same 2-3 move with several relabellings of t.
        for iso in _generateIsomorphisms( t.size(), maxIsos ):
            r = iso(t)
            source = f.embedding(0).tetrahedron().index()
            fnum = f.embedding(0).face()
            tetImage = iso.simpImage(source)
            faceImage = iso.facetPerm(source)[fnum]
            try:
                _test23single(
                        r.tetrahedron(tetImage).triangle(faceImage),
                        pach )
            except AssertionError as ae:
                print(iso)
                raise ae

    # All done!
    return count


def _test23single( face, expectedTri ):
    """
    Perform a 2-3 move on the given face, and check (among other things) that
    the result is isomorphic to expectedTri.
    """
    origTri = face.triangulation()
    tri23 = Triangulation3(origTri)
    renum = twoThree( tri23.triangle( face.index() ) )
    if not expectedTri.isIsomorphicTo(tri23):
        print(expectedTri)
        print(tri23)
        msg = "Face {}: Not isomorphic!"
        raise AssertionError( msg.format( face.index() ) )
    if ( origTri.isOriented() ) and ( not tri23.isOriented() ):
        print(tri23)
        msg = "Face {}: Failed to preserve orientation!"
        raise AssertionError( msg.format( face.index() ) )

    # Also make sure that the inverse 3-2 move gives the correct
    # triangulation, up to isomorphism.
    inv = Triangulation3(tri23)
    removedEdgeIndex = renum.underlyingEdgeIndex(-1)
    innum = threeTwo( inv.edge(removedEdgeIndex) )
    if not origTri.isIsomorphicTo(inv):
        # This test is subsumed by the more detailed isomorphisms tests below,
        # but we keep it anyway.
        print(inv)
        print( inv.detail() )
        msg = "{}: Inverse not isomorphic!"
        raise AssertionError( msg.format( face.index() ) )
    if ( tri23.isOriented() ) and ( not inv.isOriented() ):
        print(inv)
        msg = "Face {}: Inverse failed to preserve orientation!"
        raise AssertionError( msg.format( face.index() ) )
    if innum[removedEdgeIndex] is not None:
        msg = "Face {}: Inverse continues to track removed edge!"
        raise AssertionError( msg.format( face.index() ) )

    # Check that the relabellings are sensible by comparing edge
    # multiplicities in origTri and inv.
    #
    # This test is subsumed by the more detailed isomorphisms tests below, but
    # we keep it anyway.
    for i in range( origTri.countEdges() ):
        mults = _multiplicities( origTri.edge(i) )
        comMults = _multiplicities(
                inv.edge(
                    innum.underlyingEdgeIndex(
                        renum.underlyingEdgeIndex(i) ) ) )
        if mults != comMults:
            print( { k: renum.underlyingEdgeIndex(k)
                    for k in renum } )
            print( { k: innum.underlyingEdgeIndex(k)
                    for k in innum } )
            msg = "Face {}: Unmatched edge multiplicities!"
            raise AssertionError( msg.format( face.index() ) )

    # Check that the relabellings are sensible.
    if not _verifyRelab( origTri, renum, tri23, innum, inv ):
        print( { k: renum.underlyingEdgeIndex(k)
                for k in renum } )
        print( { k: innum.underlyingEdgeIndex(k)
                for k in innum } )
        msg = "Face {}: Relabellings failed!"
        raise AssertionError( msg.format( face.index() ) )

    # All done!
    return


def testTwoZero():
    total = 0
    for testSig in [ "gLLPQcdefeffpvauppb",     # 134 0-2 moves
                    "gLLPQceeffefiiaealx",      # 138 0-2 moves
                    "gvLQQcdeffeffffaafa",      # 134 0-2 moves
                    "gLLAQcdcdfffpvbbbvo" ]:    # 134 0-2 moves
        total += _test20all(testSig)
    return total


def _test20all(testSig):
    """
    Test inverse 2-0 moves on the triangulations obtained by performing 0-2
    moves on the given iso sig.
    """
    print( "0-2 and 2-0 moves on \"{}\"".format(testSig) )
    stdout.flush()
    t = Triangulation3.fromIsoSig(testSig)
    t.orient()

    count = 0
    for e in t.edges():
        deg = e.degree()
        for i in range(deg):
            for ii in range( i, deg ):
                if _test20single( t, e.index(), i, ii ):
                    count += 1

    # All done!
    return count


def _test20single( tri, edgeIndex, i, ii ):
    """
    Checks that the inverse 2-0 move correctly inverts the given 0-2 move.
    """
    tri02 = Triangulation3(tri)
    renum = _zeroTwo( tri02, edgeIndex, i, ii )
    if renum is None:
        # No 0-2 move here.
        return False
    msg = "Inverse 2-0 move {}, {}, {}: ".format( edgeIndex, i, ii )

    # Test that the inverse 2-0 move, performed on tri02 using twoZero(),
    # brings us back to a triangulation isomorphic to tri.
    inv = Triangulation3(tri02)
    innum = twoZero( inv.edge( renum.underlyingEdgeIndex(-1) ) )
    if not inv.isIsomorphicTo(tri):
        print(tri)
        print(tri02)
        print(inv)
        raise AssertionError( msg + "Not isomorphic!" )
    if ( tri02.isOriented() ) and ( not inv.isOriented() ):
        print(tri02)
        print(inv)
        raise AssertionError( msg + "Failed to preserve orientation!" )

    # Check that the relabellings are sensible.
    if not _verifyRelab( tri, renum, tri02, innum, inv ):
        print( { k: renum.underlyingEdgeIndex(k)
                for k in renum } )
        print( { k: innum.underlyingEdgeIndex(k)
                for k in innum } )
        raise AssertionError( msg + "Relabellings failed!" )

    # All done!
    return True


def _zeroTwo( tri, edgeIndex, i, ii ):
    """
    Performs a 0-2 move and returns the corresponding edge relabelling, or
    None if the requested 0-2 move is not legal.
    """
    #NOTE This implementation assumes that the two new tetrahedra introduced
    #       move02() are located at the last 2 indices.

    # Tracking edge embeddings through the requested 0-2 move is easy because
    # every existing tetrahedron survives the move.
    relab = EdgeLabelling(tri)
    if not tri.move02( tri.edge(edgeIndex), i, ii ):
        return None

    # Find the new edge that the 0-2 move just created.
    tet = tri.tetrahedron( tri.size() - 1 )
    for edgeNum in range(6):
        newEdge = tet.edge(edgeNum)
        if newEdge.degree() != 2:
            continue

        # Because the last two tetrahedra of tri were introduced by the 0-2
        # move that we just performed, the newly-introduced edge must be the
        # unique degree-2 edge that is incident to both of the last two
        # tetrahedra (this is easy to check).
        incidentTetInds = { newEdge.front().tetrahedron().index(),
                           newEdge.back().tetrahedron().index() }
        if incidentTetInds == { tri.size() - 1, tri.size() - 2 }:
            relab[-1] = newEdge.front()
            return relab
    raise AssertionError( "_zeroTwo() should never reach this point." )


def testFourFour():
    total = 0
    for testSig in [ "gLLPQcdefeffpvauppb",         # 2 4-4 moves
                    "gLLPQceeffefiiaealx",          # 2 4-4 moves
                    "gvLQQcdeffeffffaafa",          # 2 4-4 moves
                    "gLLAQcdcdfffpvbbbvo",          # 2 4-4 moves
                    "hLLzQkcdefgfggaraaavvv",       # 5 4-4 moves
                    "hLLzQkcdefgfggasaaasvs",       # 5 4-4 moves
                    "hLLzQkcdefgfggasaaavvv",       # 5 4-4 moves
                    "hLvAQkbeffgggflalaatwf",       # 5 4-4 moves
                    "ivLAPQcdefeghghhbbpbuabbv" ]:  # 5 4-4 moves
        total += _test44all(testSig)
    return total


def _test44all(testSig):
    print( "4-4 moves on \"{}\"".format(testSig) )
    stdout.flush()
    t = Triangulation3.fromIsoSig(testSig)
    t.orient()

    count = 0
    for e in t.edges():
        if _test44single(e):
            count += 1
    return count


def _test44single(edge):
    """
    Tests all possible 4-4 moves on the given edge.
    """
    t = edge.triangulation()
    #NOTE Triangulation3.has44( e, ax ) was introduced in Regina 7.4. In older
    #       versions of Regina, equivalent functionality (checking eligibility
    #       of the move, but not performing it) was provided by
    #           Triangulation3.fourFourMove( e, ax, True, False ).
    if not t.has44( edge, 0 ):
        # A 4-4 move with newAxis == 0 is available if and only if a
        # 4-4 move with newAxis == 1 is available.
        return False

    # The input edge forms one of three possible axes for an octahedron built
    # from four tetrahedra. At each such axis, we perform both possible 4-4
    # moves, and check that the isomorphism types of the resulting
    # triangulations all match up.
    isoSigSet = { i: set() for i in range(3) }
    isoSigSet[2].add( t.isoSig() )
    for newAxis in range(2):
        #NOTE Triangulation3.with44( e, ax ) was introduced in Regina 7.4.
        #       Older versions of Regina did not provide a routine to
        #       construct a new triangulation via a 4-4 move; instead, this
        #       behaviour was achieved by first constructing a copy of the
        #       triangulation, and then performing the appropriate 4-4 move on
        #       the copy.
        reg44 = t.with44( edge, newAxis )
        new44 = Triangulation3(t)
        relab = fourFour( new44.edge( edge.index() ), newAxis )

        # Test that fourFour gives the right isomorphism type, and that it
        # preserves orientation.
        move = "4-4 move {}/{}: ".format( edge.index(), newAxis )
        if not reg44.isIsomorphicTo(new44):
            print(t)
            print(reg44)
            print(new44)
            raise AssertionError( move + " Not isomorphic!" )
        if ( t.isOriented() ) and ( not new44.isOriented() ):
            raise AssertionError(
                    move + "Failed to preserve orientation!" )
        newSig = new44.isoSig()
        isoSigSet[2].add(newSig)
        isoSigSet[newAxis].add(newSig)

        # Sanity checks on the relabelling.
        if relab[ edge.index() ] is not None:
            raise AssertionError(
                    move +
                    " Relabelling continues tracking removed edge!" )
        if relab[-1] is None:
            raise AssertionError(
                    move +
                    " Relabelling fails to track new edge!" )

        # Finish filling in isoSigSet so that we can test the 4-4 moves on the
        # new edge that we just created.
        for invAxis in range(2):
            inv44 = Triangulation3(new44)
            fourFour( inv44.edge( relab.underlyingEdgeIndex(-1) ),
                     invAxis )
            isoSigSet[newAxis].add( inv44.isoSig() )
    for newAxis in range(2):
        if isoSigSet[2] != isoSigSet[newAxis]:
            raise AssertionError(
                    "Edge {}: Error with newAxis {}.".format(
                        edge.index(), newAxis ) )
    return True


def testTwoOne():
    total = 0
    for testSig in [ "dLQabccbcbv",         # 6 4-4 moves
                    "dLQabccbcsv",          # 6 4-4 moves
                    "fLAPcaccdeebgbgcv",    # 6 4-4 moves
                    "eLMkabcddbcodo",       # 8 4-4 moves
                    "eLMkabcddbcohg",       # 8 4-4 moves
                    "eLMkabcddbcoto",       # 8 4-4 moves
                    "eLMkabcddbcvag" ]:     # 8 4-4 moves
        total += _test21all(testSig)
    return total


def _test21all(testSig):
    print( "2-1 moves on \"{}\"".format(testSig) )
    stdout.flush()
    t = Triangulation3.fromIsoSig(testSig)
    t.orient()

    count = 0
    for e in t.edges():
        for edgeEnd in range(2):
            #NOTE Triangulation3.has21( e, ed ) was introduced in Regina 7.4.
            #       In older versions of Regina, equivalent functionality
            #       (checking eligibility of the move, but not performing it)
            #       was provided by
            #           Triangulation3.twoOneMove( e, ed, True, False ).
            if not t.has21( e, edgeEnd ):
                continue
            count += 1

            # Test both with the default reference labelling, and with some
            # custom reference labellings.
            _test21single( "Default", e, edgeEnd )
            trackOnlyRemovedEdge = EdgeLabelling(
                    t, { -1: e.front() } )
            _test21single( "Only removed", e, edgeEnd,
                         ( -1, trackOnlyRemovedEdge ) )
            untrackRemovedEdge = EdgeLabelling(t)
            untrackRemovedEdge.untrack( e.index() )
            _test21single( "Untrack removed", e, edgeEnd,
                         ( e.index(), untrackRemovedEdge ) )
    return count


def _test21single( name, edge, edgeEnd, data=None ):
    # Test that twoOne gives the right isomorphism type, and that it preserves
    # orientation.
    t = edge.triangulation()
    if data is None:
        removeEdgeInd = edge.index()
        edgeLab = None
        trackedIndices = [ i for i in range( t.countEdges() ) ]
        actual = Triangulation3(t)
    else:
        removeEdgeInd, edgeLab = data
        trackedIndices = edgeLab.trackedIndices()
        # Need to fully clone custom EdgeLabelling.
        edgeLab = edgeLab.clone()
        actual = edgeLab.triangulation()
    expect = t.with21( edge, edgeEnd )
    relab = twoOne( actual.edge( edge.index() ), edgeEnd, edgeLab )
    move = "{}. 2-1 move {}/{}: ".format(
            name, edge.index(), edgeEnd )
    if not actual.isIsomorphicTo(expect):
        raise AssertionError( move + "Not isomorphic!" )
    if ( t.isOriented() ) and ( not actual.isOriented() ):
        raise AssertionError(
                move + "Failed to preserve orientation!" )

    # Sanity checks on the relabelling.
    if relab[removeEdgeInd] is not None:
        raise AssertionError(
                move +
                " Relabelling continues tracking removed edge!" )
    newEdgeInd = -1 + min( 0, *trackedIndices )
    if relab[newEdgeInd] is None:
        raise AssertionError(
                move +
                " Relabelling fails to track new edge!" )

    # The new edge created by the 2-1 move should be the unique edge of degree
    # one that is incident to the new tetrahedron.
    #
    #NOTE This test assumes that the new tetrahedron is indexed last.
    newEdge = actual.edge( relab.underlyingEdgeIndex(newEdgeInd) )
    degOneEdge = None
    newTet = actual.tetrahedron( actual.size() - 1 )
    for eNum in range(6):
        candidateEdge = newTet.edge(eNum)
        if candidateEdge.degree() != 1:
            continue
        if degOneEdge is None:
            degOneEdge = candidateEdge
        else:
            raise AssertionError(
                    move +
                    " Too many degree one edges!" )
    if degOneEdge is None:
        raise AssertionError(
                move +
                " Missing degree one edge!" )
    if newEdge != degOneEdge:
        raise AssertionError(
                move +
                " New edge tracked incorrectly!" )

    # All done!
    return


def testTriGraph():
    print( "Perform moves on \"cMcabbgqs\" up to height 3." )
    stdout.flush()
    initSig = "cMcabbgqs"
    stack = [ initSig ]
    seen = { initSig }
    initSize = 2
    maxSize = 5
    counts = { "2-1": 0, "2-0": 0, "3-2": 0, "4-4": 0, "2-3": 0 }
    while stack:
        sig = stack.pop()
        tri = Triangulation3.fromIsoSig(sig)

        # Moves on edges.
        for e in tri.edges():
            if e.degree() == 1:
                for edgeEnd in {0,1}:
                    newTri = Triangulation3(tri)
                    relab = twoOne( newTri.edge( e.index() ), edgeEnd )
                    if relab is None:
                        continue
                    counts["2-1"] += 1

                    # If we haven't seen the new triangulation before (up to
                    # combinatorial isomorphism), then add it to the stack.
                    newSig = newTri.isoSig()
                    if newSig not in seen:
                        seen.add(newSig)
                        stack.append(newSig)
            elif e.degree() == 2:
                newTri = Triangulation3(tri)
                relab = twoZero( newTri.edge( e.index() ) )
                if relab is None:
                    continue
                counts["2-0"] += 1

                # If we haven't seen the new triangulation before (up to
                # combinatorial isomorphism), then add it to the stack.
                newSig = newTri.isoSig()
                if newSig not in seen:
                    seen.add(newSig)
                    stack.append(newSig)
            elif e.degree() == 3:
                newTri = Triangulation3(tri)
                relab = threeTwo( newTri.edge( e.index() ) )
                if relab is None:
                    continue
                counts["3-2"] += 1

                # If we haven't seen the new triangulation before (up to
                # combinatorial isomorphism), then add it to the stack.
                newSig = newTri.isoSig()
                if newSig not in seen:
                    seen.add(newSig)
                    stack.append(newSig)
            elif e.degree() == 4:
                for newAxis in {0,1}:
                    newTri = Triangulation3(tri)
                    relab = fourFour( newTri.edge( e.index() ), newAxis )
                    if relab is None:
                        continue
                    counts["4-4"] += 1

                    # If we haven't seen the new triangulation before (up to
                    # combinatorial isomorphism), then add it to the stack.
                    newSig = newTri.isoSig()
                    if newSig not in seen:
                        seen.add(newSig)
                        stack.append(newSig)

        # 2-3 moves (only if such moves would not exceed the given height).
        if tri.size() == maxSize:
            continue
        for f in tri.triangles():
            newTri = Triangulation3(tri)
            relab = twoThree( newTri.triangle( f.index() ) )
            if relab is None:
                continue
            counts["2-3"] += 1

            # If we haven't seen the new triangulation before (up to
            # combinatorial isomorphism), then add it to the stack.
            newSig = newTri.isoSig()
            if newSig not in seen:
                seen.add(newSig)
                stack.append(newSig)

    # Print move counts.
    print()
    print( "Tested: {} {}; {} {}; {} {}; {} {}; {} {}.".format(
        counts["2-1"], "2-1 moves",
        counts["2-0"], "2-0 moves",
        counts["3-2"], "3-2 moves",
        counts["4-4"], "4-4 moves",
        counts["2-3"], "2-3 moves" ) )
    return sum( counts.values() )


if __name__ == "__main__":
    RandomEngine.reseedWithHardware()
    availableTests = [ "23", "20", "44", "21", "graph" ]
    testNames = parseTestNames( argv[1:], availableTests )

    # Test 2-3 and 3-2 moves.
    if "23" in testNames:
        longName = "2-3 and 3-2 moves"
        runNamedTestSuite( longName, testTwoThree )

    # Test 2-0 moves.
    if "20" in testNames:
        longName = "2-0 moves"
        runNamedTestSuite( longName, testTwoZero )

    # Test 4-4 moves.
    if "44" in testNames:
        longName = "4-4 moves"
        runNamedTestSuite( longName, testFourFour )

    # Test 2-1 moves.
    if "21" in testNames:
        longName = "2-1 moves"
        runNamedTestSuite( longName, testTwoOne )

    # Perform a bunch of moves to check that we don't get any exceptions.
    if "graph" in testNames:
        longName = "Triangulation graph"
        runNamedTestSuite( longName, testTriGraph )

    # If we make it here, then all tests passed.
    allTestsPassedMessage(testNames)
