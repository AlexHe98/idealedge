"""
Construct triangulations of orientable Seifert fibre spaces.
"""
from sys import argv
from regina import *
from surftri import surface
from test import parseTestNames, doTest, allTestsPassedMessage


class OrientableBundle:
    """
    An oriented triangulation of an (orientable) bundle over a surface.

    The bundle can either be a circle bundle or an I-bundle.

    The 3-manifold triangulation T of the bundle is constructed from a
    2-manifold triangulation S by assigning a triangular prism (see the
    TriPrism class) to each triangle of S, and gluing the prisms together
    accordingly. For triangle i of S, the associated triangular prism can
    be accessed using OrientableBundle.triPrism(i); moreover, edge e of
    triangle i corresponds to boundary square number e of the associated
    prism.
    """
    def __init__( self, surf, circleBundle ):
        """
        Constructs an oriented triangulation of an (orientable) bundle over
        the given base surface.

        If circleBundle is True, then this will be a circle bundle.
        Otherwise, it will be an I-bundle.
        """
        self._base = surf
        self._circleBundle = circleBundle

        # Create triangular prisms that give bundles over each triangle in base,
        # and translate the gluings of base into gluings of these prisms.
        self._tri = Triangulation3()
        self._prism = [ TriPrism( self._tri, circleBundle )
                       for _ in range( surf.size() ) ]
        for edge in surf.edges():
            if edge.isBoundary():
                continue

            # Get gluing data from the surface triangulation.
            myFaceIndex = edge.front().triangle().index()
            myEdgeNum = edge.front().edge()
            yourFaceIndex = edge.back().triangle().index()
            gluing = edge.front().triangle().adjacentGluing(myEdgeNum)

            # Translate this into a gluing of prisms.
            myPrism = self._prism[myFaceIndex]
            yourPrism = self._prism[yourFaceIndex]
            myPrism.orientedJoin( myEdgeNum, yourPrism, gluing )

        # All done!
        return

    def triangulation(self):
        """
        Returns the underlying 3-manifold triangulation of this bundle.
        """
        return self._tri

    def base(self):
        """
        Returns the triangulation of the base surface of this bundle.
        """
        return self._base

    def triPrism( self, i ):
        """
        Returns the triangular prism corresponding to triangle i of the base
        surface.
        """
        return self._prism[i]

    def isCircleBundle(self):
        """
        Is this a circle bundle?
        """
        return self._circleBundle

    def isIBundle(self):
        """
        Is this an I-bundle?
        """
        return not self._circleBundle

    def toggle(self):
        """
        Toggles whether this is a circle bundle or an I-bundle.
        """
        for prism in self._prism:
            prism.toggleSolidTorus()
        self._circleBundle = not self._circleBundle
        return

    # End of OrientableBundle class.
    pass


class TriPrism:
    r"""
    An oriented triangular prism within a triangulation.

    The vertical boundary of such a triangular prism consists of three
    squares. For s in {0,1,2}, square number s is the one that is opposite
    tetrahedron s of the triangulation.

    Each boundary square is triangulated by a pair of triangles. The choice
    of direction for the diagonal edge separating the two triangles is
    specified by a slope: either +1 or -1. Within each square, the triangles
    are numbered either 0 (left triangle in the diagrams below) or 1 (right
    triangle in the diagrams below).

                        Slope +1        Slope -1
                         •----•          •----•
                         |0 2/|          |\2 0|
                         |  /1|          |1\  |
                         |1/  |          |  \1|
                         |/2 0|          |0 2\|
                         •----•          •----•

    Details about boundary square number s can via the following three pieces
    of data:
    --> squareSlope(s) gives the slope of the diagonal of square s.
    --> squareTet( s, t ) gives the tetrahedron that provides triangle t of
        square s.
    --> squareRoles( s, t ) gives a permutation that specifies the roles
        played by the vertices of squareTet( s, t ) with respect to the
        boundary square. In detail:
        --- For v in {0,1,2}, squareRoles( s, t )[v] gives the vertex
            labelled v in the above diagrams.
        --- squareRoles( s, t )[3] gives the vertex that is opposite the
            triangle.
    For each s in {0,1,2} and t in {0,1}, it is guaranteed that
        squareSlope(s) * squareRoles( s, t ).sign() == -1.

    As shown in the above diagram, for triangle t of square s, the role of
    each edge is determined by its endpoints:
    --> If the endpoints play roles 0 and 1, then the edge is vertical.
    --> If the endpoints play roles 0 and 2, then the edge is horizontal.
    --> If the endpoints play roles 1 and 2, then the edge is diagonal.
    We number the vertical edges by 0, 1 or 2 depending on which tetrahedron
    they are incident to. Thus, for each i in {0,1,2}, vertical edge i is
    opposite boundary square i.

    This class provides the flipSlope(s) method, which flips the slope of
    boundary square s by either adding a tetrahedron (via layering) or
    removing a tetrahedron (if the slope was already previously flipped via
    layering). The outputs for squareTet() and squareRoles() are
    automatically adjusted to account for any such layerings.
    """
    def __init__( self, tri, solidTorus=True ):
        """
        Constructs a triangular prism and inserts it into the given
        3-manifold triangulation tri.

        If solidTorus is True (the default), then the top and bottom will be
        glued together so that the triangular prism forms a solid torus.
        Otherwise, the top and bottom will be left unglued and the triangular
        prism will form a 3-ball.
        """
        self._tri = tri

        # Build the triangular prism from three "core" tetrahedra (as opposed
        # to tetrahedra that might be introduced later via layering).
        self._coreTet = tri.newTetrahedra(3)
        self._coreTet[0].join( 1, self._coreTet[1], Perm4(2,3) )
        self._coreTet[1].join( 3, self._coreTet[2], Perm4(2,3) )
        if solidTorus:
            self._glueTopBot()

        # Build the boundary square data.
        #NOTE We can manually check that the chosen values below satisfy the
        #   requirement that slope * sign == -1.
        self._squareOrigSlope = [ 1, -1, 1 ]
        self._squareOrigTetIndex = [
                [ 2, 1 ],
                [ 0, 2 ],
                [ 1, 0 ] ]
        self._squareOrigRoles = [
                [ Perm4(2,3,1,0), Perm4(3,1,2,0) ],
                [ Perm4(1,0,3,2), Perm4(2,3,0,1) ],
                [ Perm4(1,3,0,2), Perm4(1,0,2,3) ] ]

        # If we introduce layerings later on, then this will change the
        # original boundary square data recorded above.
        self._layerTet = [ None, None, None ]
        self._layerRoles = [ None, None, None ]

        # We will need to keep track of which boundary squares of this prism
        # are glued to other prisms.
        self._isGlued = [ False, False, False ]
        return

    def _glueTopBot(self):
        self._coreTet[2].join( 3, self._coreTet[0], Perm4(1,2,3,0) )

    def triangulation(self):
        """
        Returns the triangulation containing this triangular prism.
        """
        return self._tri

    def isSolidTorus(self):
        """
        Are the top and bottom glued, so that this triangular prism forms a
        solid torus?
        """
        adj = self._coreTet[2].adjacentTetrahedron(3)
        return ( adj is not None )

    def toggleSolidTorus(self):
        """
        Toggles whether this triangular prism forms a solid torus by either
        gluing or ungluing the top and bottom.

        This routine returns the new status: True if this triangular prism
        becomes a solid torus after the toggle, and False otherwise.
        """
        if self.isSolidTorus():
            self._coreTet[2].unjoin(3)
            return False
        else:
            self._glueTopBot()
            return True
        raise AssertionError( "This should be unreachable." )

    def squareSlope( self, s ):
        """
        Returns the current slope of the diagonal edge of boundary square s
        of this triangular prism.
        """
        origSlope = self._squareOrigSlope[s]
        if self._layerTet[s] is None:
            return origSlope
        return -origSlope

    def squareTet( self, s, t ):
        """
        Returns the tetrahedron that provides triangle t of boundary square s
        of this triangular prism.
        """
        if self._layerTet[s] is None:
            return self._coreTet[ self._squareOrigTetIndex[s][t] ]
        return self._layerTet[s]

    def squareRoles( self, s, t ):
        """
        Returns the permutation describing the roles of the vertices of
        self.squareTet(s,t).

        See the class notes for details.
        """
        if self._layerTet[s] is None:
            return self._squareOrigRoles[s][t]
        return self._layerRoles[s][t]

    def isSquareGlued( self, s ):
        """
        Is boundary square s of this prism glued to another prism?
        """
        return self._isGlued[s]

    def flipSlope( self, square ):
        """
        Flips the sign of the slope of the given square's diagonal edge.

        This routine returns the new slope that results from performing the
        flip.

        Pre-condition:
        --> self.isSquareGlued(square) is False.
        """
        doomed = self._layerTet[square]
        if doomed is not None:
            # This is the easy case: we can flip the slope by simply removing
            # the previously layered tetrahedron.
            self._layerTet[square] = None
            self._layerRoles[square] = None
            self._tri.removeTetrahedron(doomed)

            # We have reverted to the original slope.
            return self._squareOrigSlope[square]

        # No previous layering, so we need to introduce such a layering.
        #
        #               Slope +1          Slope -1
        #             Roles sign -1     Roles sign +1
        #                •----•            •----•
        #                |0 2/|            |\2 0|
        #                |  /1|            |1\  |
        #                |1/  |            |  \1|
        #                |/2 0|            |0 2\|
        #                •----•            •----•
        #
        # We will flip the diagonal by layering a tetrahedron. The following
        # vertex labelling ensures that the layering preserves orientation:
        #
        #               2•----•1          0•----•2
        #                |\   |            |   /|
        #                | \  |            |  / |
        #                |  \ |            | /  |
        #                |   \|            |/   |
        #               0•----•3          3•----•1
        #
        if self._squareOrigSlope[square] == 1:
            # Sanity check: roles sign is already -1, so the signs of our
            # gluing permutations should be +1 (as they indeed are).
            relativeGluing = [ Perm4(2,0,1,3), Perm4(3,1,0,2) ]

            # Slope will become -1, so the sign of the new roles permutations
            # needs to be +1 (as they indeed are).
            self._layerRoles[square] = [ Perm4(0,2,3,1), Perm4(1,3,2,0) ]
        else:
            # Sanity check: roles sign is +1, so the signs of our gluing
            # permutations should be -1 (as they indeed are).
            relativeGluing = [ Perm4(3,0,1,2), Perm4(2,1,0,3) ]

            # Slope will become +1, so the sign of the new roles permutations
            # needs to be -1 (as they indeed are).
            self._layerRoles[square] = [ Perm4(0,3,2,1), Perm4(1,2,3,0) ]

        # OK, actually perform the layering.
        layerTet = self._tri.newTetrahedron()
        self._layerTet[square] = layerTet
        for t in range(2):
            coreTet = self._coreTet[ self._squareOrigTetIndex[square][t] ]
            origRoles = self._squareOrigRoles[square][t]
            coreFace = origRoles[3]
            gluing = relativeGluing[t] * origRoles.inverse()
            coreTet.join( coreFace, layerTet, gluing )
        return -self._squareOrigSlope[square]

    def orientedJoin( self, s, otherPrism, gluing ):
        """
        Joins boundary square s of this prism to some boundary square of
        another prism in the same triangulation.

        You may join a boundary square of this prism to some different
        boundary square of the same prism (i.e., otherPrism is allowed to be
        this prism), though you cannot join a boundary square to itself.

        If the triangulation is currently oriented, then the requested gluing
        will preserve this orientation.

        The rationale for this routine is that if we associate a collection
        of prisms to the triangles of a 2-manifold triangulation S, then
        joining all the prisms together using the same gluing permutations as
        in S will produce an oriented 3-manifold triangulation T such that:
        --> if isSolidTorus() is True for every prism, then T will be the
            orientable circle bundle over S; and
        --> if isSolidTorus() is False for every prism, then T will be the
            orientable I-bundle over S.

        Pre-condition:
        --> This prism and otherPrism belong to the same triangulation.
        --> Boundary square s of this prism is not currently glued to
            anything.
        --> The corresponding boundary square of otherPrism (i.e., boundary
            square gluing[s] of otherPrism) is likewise not currently glued
            to anything.
        --> We are not attempting to glue a boundary square to itself (i.e.,
            we do not have both self == otherPrism and gluing[s] == s).

        Parameters:
        --> s           The boundary square of this prism that will be glued
                        to otherPrism; s must be in {0,1,2}.
        --> otherPrism  The other prism that will be glued to boundary square
                        s of this prism.
        --> gluing      A permutation that describes how the vertical edges
                        of this prism will map to vertical edges otherPrism
                        across the new gluing.
        """
        sOther = gluing[s]

        # Check some of the preconditions.
        if self.isSquareGlued(s) or otherPrism.isSquareGlued(sOther):
            raise ValueError(
                    "Can only glue squares that are currently unglued." )

        # We need the two squares on either side of the gluing to have
        # opposite slopes.
        if self.squareSlope(s) == otherPrism.squareSlope(sOther):
            # Flip slope by removing tetrahedra whenever possible.
            if ( ( self._layerTet[s] is not None ) or
                ( otherPrism._layerTet[sOther] is None ) ):
                self.flipSlope(s)
            else:
                otherPrism.flipSlope(sOther)

        # Recall that the boundary squares are labelled as follows:
        #
        #               Slope +1          Slope -1
        #             Roles sign -1     Roles sign +1
        #                •----•            •----•
        #                |0 2/|            |\2 0|
        #                |  /1|            |1\  |
        #                |1/  |            |  \1|
        #                |/2 0|            |0 2\|
        #                •----•            •----•
        #
        for t in range(2):
            if gluing.sign() == 1:
                # Gluing swaps top and bottom, not left and right.
                tOther = t
            else:
                # Gluing swaps left and right.
                tOther = 1-t
            roles = self.squareRoles( s, t )
            rolesOther = otherPrism.squareRoles( sOther, tOther )

            # Perform the gluing on triangle t of square s of this prism.
            # Opposite slopes implies that roles and rolesOther have opposite
            # signs, so our choice of tetGluing permutation does indeed have
            # sign -1.
            myTet = self.squareTet( s, t )
            myFace = roles[3]
            yourTet = otherPrism.squareTet( sOther, tOther )
            tetGluing = rolesOther * roles.inverse()
            myTet.join( myFace, yourTet, tetGluing )

        # All done!
        self._isGlued[s] = True
        otherPrism._isGlued[sOther] = True
        return

    # End of TriPrism class.
    pass


if __name__ == "__main__":
    availableTests = [ "bundle",
                      "prism" ]
    testNames = parseTestNames( argv[1:], availableTests )

    # Test OrientableBundle class.
    if "bundle" in testNames:
        print( "+------------------------+" )
        print( "| OrientableBundle class |" )
        print( "+------------------------+" )

        # Check that we get oriented triangulations of the correct
        # 3-manifolds.
        #NOTE Regina is not always consistent with how it chooses from the
        #   many alternative names that these manifolds have.
        expectedNames = [
                "KB/n2 x~ S1",
                "Non-or, g=2 + 1 puncture/n2 x~ S1",
                "Non-or, g=2 + 2 punctures/n2 x~ S1",
                "RP3 # RP3",
                "M/n2 x~ S1",
                "Non-or, g=1 + 2 punctures/n2 x~ S1",
                "S2 x S1",
                "SFS [D: (1,1)]",
                "SFS [A: (1,1)]",
                "T x S1",
                "Or, g=1 + 1 puncture x S1",
                "Or, g=1 + 2 punctures x S1",
                "Or, g=2 x S1",
                "Or, g=2 + 1 puncture x S1",
                "Or, g=2 + 2 punctures x S1" ]
        nameIndex = -1
        for genus in range( -2, 3 ):
            for boundaries in range(3):
                print( "g={}, b={}.".format( genus, boundaries ) )
                nameIndex += 1
                expectedName = expectedNames[nameIndex]
                base = surface( genus, boundaries )

                # Circle bundle.
                bundle = OrientableBundle( base, True )
                tri = bundle.triangulation()
                doTest( "Oriented?", True, tri.isOriented() )
                actual = StandardTriangulation.recognise(tri).manifold()
                doTest( "Manifold?", expectedName, actual.name() )
                for edge in base.edges():
                    if not edge.isBoundary():
                        continue
                    faceIndex = edge.front().triangle().index()
                    prism = bundle.triPrism(faceIndex)
                    square = edge.front().edge()
                    doTest(
                            "Square {}({}) glued?".format(
                                faceIndex, square ),
                            False, prism.isSquareGlued(square) )
                print()

    # Test TriPrism class.
    if "prism" in testNames:
        print( "+----------------+" )
        print( "| TriPrism class |" )
        print( "+----------------+" )

        def _testTri( tri, isSolidTorus ):
            """
            Test basic properties of the triangulation.
            """
            doTest( "Solid torus?", isSolidTorus, tri.isSolidTorus() )
            doTest( "Ball?", not isSolidTorus, tri.isBall() )
            doTest( "Oriented?", True, tri.isOriented() )
            return

        def _testSlopeSign( prism, s, tri, isSolidTorus ):
            """
            Test that the slope-sign constraint holds both before and after
            flipping square s of the given prism.
            """
            # Before flipping.
            oldSlope = prism.squareSlope(s)
            print( "Old slope of square {}: {}.".format(
                s, oldSlope ) )
            for t in range(2):
                doTest( "Slope-sign constraint for triangle {}.".format(t),
                       -1, oldSlope * prism.squareRoles(s,t).sign() )

            # After flipping.
            newSlope = prism.flipSlope(s)
            doTest( "Flip slope of square {}.".format(s),
                    -oldSlope, newSlope )
            doTest( "New slope of square {}.".format(s),
                    -oldSlope, prism.squareSlope(s) )
            _testTri( tri, isSolidTorus )
            for t in range(2):
                doTest( "Slope-sign constraint for triangle {}.".format(t),
                        -1, newSlope * prism.squareRoles(s,t).sign() )
            print()
            return

        # Test both solid-torus and non-solid-torus constructions.
        for isSolidTorus in ( True, False ):
            tri = Triangulation3()
            prism = TriPrism( tri, isSolidTorus )
            print( "------------------------------------------------" )
            _testTri( tri, isSolidTorus )
            print( "------------------------------------------------" )
            print()

            # Flip everything one by one, and then unflip, and test that the
            # slope-sign constraint is preserved all the way through.
            for _ in range(2):
                for s in range(3):
                    _testSlopeSign( prism, s, tri, isSolidTorus )

    # If we make it here, then all tests passed.
    allTestsPassedMessage(testNames)
