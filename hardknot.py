"""
Try to generate a hard diagram of a composite knot by randomly composing
prime knots using the overlaying construction.
"""
from sys import argv, stdout
from timeit import default_timer
from multiprocessing import Pool, TimeoutError
import snappy
from regina import *
from hardknotimpl import attemptHardComposite


def randomHardComposite( numSummands, workers, verbose=True ):
    """
    Randomly generate a hard diagram of a composite knot by repeatedly
    composing n random prime knots, where n = numSummands.

    Since generating the desired hard diagram often requires many attempts
    and a significant amount of time, it is often helpful to run this routine
    with multiple concurrent worker processes. It may also be helpful to
    regularly print progress updates using verbose mode (which is switched on
    by default).
    """
    start = default_timer()
    if verbose:
        prev = start
    total = 0
    processed = 0
    unfinished = dict()
    maxSize = 2*workers
    msg = "Time: {:.6f}. Generated: {}. Processed: {}."
    with Pool( processes=workers ) as pool:
        # Generate random composite knots, overlay, and check whether this
        # leads to a hard diagram.
        for _ in range(maxSize):
            unfinished[total] = pool.apply_async(
                    attemptHardComposite, args = (numSummands,) )
            total += 1
        while unfinished:
            if len(unfinished) < maxSize:
                unfinished[total] = pool.apply_async(
                        attemptHardComposite, args = (numSummands,) )
                total += 1

            # Check for new results.
            nowFinished = set()
            for i in unfinished:
                try:
                    output = unfinished[i].get( timeout=0.001 )
                except TimeoutError:
                    pass
                else:
                    # Got a new result!
                    nowFinished.add(i)
                    processed += 1

                    # Is it the hard diagram we want?
                    if output is not None:
                        names, compPD = output
                        comp = snappy.Link(compPD)
                        print( msg.format(
                            default_timer() - start, total, processed ) )
                        return names, compPD, comp
            for i in nowFinished:
                del unfinished[i]

            # In verbose mode, print periodic updates.
            if verbose:
                time = default_timer()
                if time - prev > 10:
                    prev = time
                    print( msg.format(
                        time - start, total, processed ) )
                    stdout.flush()

    # Should never reach this point.
    raise RuntimeError()


if __name__ == "__main__":
    numSummands = int( argv[1] )
    workers = int( argv[2] )
    names, pd, comp = randomHardComposite( numSummands, workers )
    for n in names:
        print(n)
    print()
    print(comp)
    print()

    # View the knot diagram using the PLink viewer.
    #NOTE This doesn't seem to work on all machines.
    comp.view().window.mainloop()

    # Print knot signature.
    print( Link.fromPD(pd).knotSig() )
