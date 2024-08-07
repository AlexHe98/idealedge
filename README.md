Ideal loops
-----------

This repository contains source code for working with *ideal loops*: loops
that are embedded in the interior of a 3-manifold triangulation as a sequence
of edges. The rationale is that a *loop-ideal triangulation* (i.e., a
triangulation endowed with an ideal loop) represents the 3-manifold with
boundary given by deleting a small regular neighbourhood of the ideal loop.

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

— *Alex He (a.he@uqconnect.edu.au)*
