Edge-ideal triangulations
-------------------------

This repository implements algorithms which rely on the notion of an
*edge-ideal triangulation*: that is, a 3-manifold triangulation *T*
containing a collection of distinguished loops, called *ideal loops*, which
are embedded as sequences of internal edges of *T*.

The algorithms implemented here are underpinned by theory which was developed
in the following joint work with *Eric Sedgwick* and *Jonathan Spreer*:
- A practical algorithm for knot factorisation.\
    [DOI:10.4230/LIPIcs.SoCG.2025.55](https://doi.org/10.4230/LIPIcs.SoCG.2025.55)\
    [arXiv:2504.03942](https://arxiv.org/abs/2504.03942)
- Practical bounded orientable Seifert fibred space recognition.\
    *In preparation.*

This source code depends on [Regina](https://regina-normal.github.io/). The
oldest compatible version of Regina is 7.4. For users who are (for whatever
reason) unable to upgrade to Regina 7.4 or later, there are some comments
throughout the code to help with porting to older versions of Regina.

An edge-ideal triangulation represents a 3-manifold *M* which has
"edge-ideal" torus boundary components given by drilling out each of the
ideal loops. One of the advantages of this idea is that a normal surface in
an edge-ideal triangulation must intersect ideal loops traversely, and must
therefore correspond to a properly embedded surface in *M* that intersects
the edge-ideal boundary tori in curves of a prescribed slope (namely, the
slope of the meridian of a solid torus neighbourhood of an ideal loop).

As suggested by the paper titles listed above, there are currently two main
applications for the techniques implemented here:
- An algorithm for decomposing a knot into its prime summands (see the
    ``decompose()`` routine in ``decomposeknot.py``). The output primes are
    always given as edge-ideal triangulations, but the input knot is more
    flexible: the input could be an edge-ideal triangulation, but could also
    be a Regina ``Link`` object. Naturally, this algorithm specialises in the
    following ways:
    + *Unknot recognition*: A knot is the unknot if and only if its prime
        factorisation is empty.
    + *Prime knot recognition*: A nontrivial knot is prime if and only if its
        prime factorisation contains exactly one knot.
    + *Composite knot recognition*: A knot is composite if and only if its
        prime factorisation contains at least two knots.
- An algorithm for recognising bounded orientable Seifert fibred spaces (see
    the ``recogniseSFS()`` routine in ``recsfs.py``). This algorithm is
    already fully functional, but further improvements are still in
    development. For knots, Seifert fibred space recognition specialises to
    torus knot recognition (see the ``recogniseTorusKnot()`` routine, also in
    ``recsfs.py``).

The main scripts in this repository are the following:
- ``decomposeknot.py``: Contains the ``decompose()`` routine.
- ``recsfs.py``: Contains the ``recogniseSFS()`` and ``recogniseTorusKnot()``
    routines.
- ``loop.py``: Implements the ``IdealLoop`` class.
- ``triloops.py``: Implements the ``EdgeIdealTriangulation`` class.
- ``embed.py``: Implements routines for converting a Regina ``Link`` object
    into an ``EdgeIdealTriangulation``.
- ``idealedge.py``: Implements routines for crushing a normal surface, while
    keeping track of not just how the ambient triangulation changes, but also
    how the ideal loops change.
- ``demo.py``: Runs a live demonstration of the ``decompose()`` routine for
    knots, either with a randomly-generated hard diagram of a composite knot,
    or with a user-provided knot signature.

An important test of the performance of our knot decomposition algorithm is on
hard diagrams of composite knots. Here, by *hard*, we mean that the diagram is
diagrammatically prime, and cannot be simplified using SnapPy's global
simplification heuristic. Our code for generating such hard diagrams, which
may be of independent interest, is available in the ``hardknot/`` directory.

This repository includes, in the ``experiments/knots/`` directory, the
following scripts for running computational experiments:
- ``experiments/knots/main.py``: Runs the ``decompose()`` routine on all
    knots from a given collection of knot tables.
- ``experiments/knots/sample.py``: Runs the ``decompose()`` routine on a
    random sample of knots from a given collection of knot tables.
- ``experiments/knots/composite.py``: Runs the ``decompose()`` routine on a
    collection of composite knots constructed by composing knots that are
    randomly sampled from a given collection of knot tables.
- ``experiments/knots/torus.py``: Runs the ``decompose()`` routine on all
    torus knots with crossing number in a given interval.
- ``experiments/knots/s3edges.py``: Runs the ``decompose()`` routine on all
    knots that appear as edges in a given collection of one-vertex
    triangulations of the 3-sphere.
- ``experiments/knots/sigs.py``: Runs the ``decompose()`` routine on all
    knots given by knot signatures in a given dataset.

One of the implementation challenges with edge-ideal triangulations is
keeping track of the ideal loops as we modify the ambient triangulation (say,
through local moves, or through crushing). In particular, our implementations
for tracking edges through local moves, available in the ``retriangulate/``
directory, may be of independent interest.

— *Alex He (a.he@uqconnect.edu.au)*, 31 Aug 2026
