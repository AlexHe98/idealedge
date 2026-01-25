Ideal loops
-----------

This repository contains source code for working with *ideal loops*: loops
that are embedded in the interior of a 3-manifold triangulation as a sequence
of edges.

The rationale is that an *edge-ideal triangulation* (i.e., a triangulation
endowed with an ideal loop) represents a 3-manifold *M* with a torus boundary
component *B* given by deleting a small regular neighbourhood *N* of the
ideal loop. One of the advantages of this notion is that a normal surface in
an edge-ideal triangulation must intersect the ideal loop transversely, and
therefore corresponds to a surface that intersects the boundary torus *B* in
curves with a prescribed slope (specifically, the curves correspond to the
meridian of the solid torus *N*).

Currently, the main application is an algorithm for decomposing a knot into
its prime summands. The output primes are always given as ideal loops, but
the input knot is more flexible: the input could be an ideal loop, but could
also be a Regina Link object. This knot decomposition algorithm was designed
in joint work with *Eric Sedgwick* and *Jonathan Spreer*.

The main scripts in this repository are the following:
- ``decomposeknot.py``: Contains the decompose() routine, which implements
    the aforementioned knot decomposition algorithm.
- ``loop.py``: Implements the IdealLoop class.
- ``embed.py``: Implements routines for converting a Regina Link object into
    an ideal loop.
- ``idealedge.py``: Contains the decomposeAlong() routine, which crushes a
    normal surface, while keeping track of not just how the triangulation
    changes, but also how the ideal loop changes.

This repository also includes, in the ``experiments/knots/`` directory, the
following scripts for running computational experiments:
- ``experiments/knots/main.py``: Runs the decompose() routine on all knots
    from a given collection of knot tables.
- ``experiments/knots/sample.py``: Runs the decompose() routine on a random
    sample of knots from a given collection of knot tables.
- ``experiments/knots/composite.py``: Runs the decompose() routine on a
    collection of composite knots constructed by composing knots that are
    randomly sampled from a given collection of knot tables.
- ``experiments/knots/torus.py``: Runs the decompose() routine on all torus
    knots with crossing number in a given interval.
- ``experiments/knots/s3edges.py``: Runs the decompose() routine on all knots
    that appear as edges in a given collection of one-vertex triangulations of
    the 3-sphere.
- ``experiments/knots/sigs.py``: Runs the decompose() routine on all knots
    given by knot signatures in a given dataset.

— *Alex He (a.he@uqconnect.edu.au)*
