"""
Generate triangulations of bounded orientable 3-manifolds other than a
Seifert fibred space or hyperbolic 3-manifold.
"""
from dataclasses import dataclass
from pathlib import Path
from timeit import default_timer
from regina import *
import snappy
from experiments.sfs.io import hardSFSPath, readSFSLine, parseSFSFileName


MAX_CUSPS = 4


def _simplifyHard(tri):
    # Try really hard to simplify.
    simplifiedNow = True
    while simplifiedNow:
        simplifiedNow = tri.simplify()
        if not simplifiedNow:
            simplifiedNow = tri.simplify()
    return


def _randomHyp(numBdries):
    mfd = snappy.OrientableCuspedCensus(num_cusps=numBdries).random()
    tri = Triangulation3( mfd.triangulation_isosig(decorated=False) )
    tri.idealToFinite()
    tri.simplify()
    tri.minimiseBoundary()
    _simplifyHard(tri)

    # Apply a random isomorphism.
    return Isomorphism3.random( tri.size() )(tri)


@dataclass
class _Gluing:
    """
    A gluing between tetrahedra.
    """
    myTeti: int
    myFace: int
    yourTeti: int
    gluingPerm: Perm4

    def performOn(self, tri):
        """
        Performs this gluing on the given triangulation.
        """
        tri.tetrahedron(self.myTeti).join( self.myFace,
                                          tri.tetrahedron(self.yourTeti),
                                          self.gluingPerm )
        return


if __name__ == "__main__":
    RandomEngine.reseedWithHardware()
    start = default_timer()

    # Only need to do this once, so just hard-code the path where we store
    # the data.
    with open( "data/nonsfs/sample.txt", "x" ) as outputFile:
        for file in Path( hardSFSPath() ).iterdir():
            if file.is_file():
                isBaseOrbl, genus, numBdries = parseSFSFileName(file.name)
                if numBdries < 2:
                    continue

                # We have a collection of hard SFS with 2 or more boundary
                # tori. Attach a hyperbolic 3-manifold, but leaving at least
                # one of the SFS boundaries exposed, so that recogniseSFS()
                # has a chance of finding some surface to crush.
                if isBaseOrbl:
                    baseClass = SFSpace.Class.bo1
                else:
                    baseClass = SFSpace.Class.bn2
                with open( hardSFSPath(
                    isBaseOrbl, genus, numBdries ) ) as sfsFile:
                    for line in sfsFile:
                        neoSig, _fibresIgnored = readSFSLine(line)
                        sfs = Triangulation3.fromSig(neoSig)
                        sfs = Isomorphism3.random( sfs.size() )(sfs)

                        # Pick the number of cusps in the hyperbolic piece to
                        # attach.
                        numCuspsToAttach = 1 + RandomEngine.rand(
                                min( MAX_CUSPS, numBdries - 1 ) )
                        numCusps = numCuspsToAttach + RandomEngine.rand(
                                1 + MAX_CUSPS - numCuspsToAttach )
                        hyp = _randomHyp(numCusps)
                        assert ( numCuspsToAttach <
                                sfs.countBoundaryComponents() )
                        assert ( numCuspsToAttach <=
                                hyp.countBoundaryComponents() )

                        # Figure out the gluings.
                        gluings = []
                        for i in range(numCuspsToAttach):
                            sfsBdry = sfs.boundaryComponent(i)
                            hypBdry = hyp.boundaryComponent(i)
                            sfsEdge = sfsBdry.edge(0)
                            hypEdge = hypBdry.edge(0)

                            # Gluing coming from the front.
                            sfsFront = sfsEdge.front()
                            hypFront = hypEdge.front()
                            sfsTeti = sfsFront.tetrahedron().index()
                            sfsFace = sfsFront.vertices()[3]
                            hypTeti = ( sfs.size() +
                                       hypFront.tetrahedron().index() )
                            perm = ( hypFront.vertices() *
                                    sfsFront.vertices().inverse() )
                            gluings.append( _Gluing(
                                sfsTeti, sfsFace, hypTeti, perm ) )

                            # Gluing coming from the back.
                            sfsBack = sfsEdge.back()
                            hypBack = hypEdge.back()
                            sfsTeti = sfsBack.tetrahedron().index()
                            sfsFace = sfsBack.vertices()[2]
                            hypTeti = ( sfs.size() +
                                       hypBack.tetrahedron().index() )
                            perm = ( hypBack.vertices() *
                                    sfsBack.vertices().inverse() )
                            gluings.append( _Gluing(
                                sfsTeti, sfsFace, hypTeti, perm ) )

                        # Perform the gluings.
                        sfs.insertTriangulation(hyp)
                        for gluing in gluings:
                            gluing.performOn(sfs)

                        # Do some sanity checks before saving this example.
                        assert sfs.isValid()
                        assert sfs.isOrientable()
                        assert ( sfs.countBoundaryComponents() ==
                                numBdries + numCusps - 2*numCuspsToAttach )
                        outputFile.write( f"{sfs.neoSig()}\n" )
    print( "Time: {:.6f}".format( default_timer() - start ) )
