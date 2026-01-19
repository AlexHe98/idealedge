"""
Helper routines for test suites.
"""


def parseTestNames( arguments, availableTests ):
    """
    Parses a list of arguments into a set of test names chosen from the list
    of available tests.

    In addition to the names explicitly listed in availableTests, this
    routine also accepts "all" as a valid test name, and interprets "all" as
    a request to perform all tests listed in availableTests.

    If the list of arguments is empty, or if there is any single argument
    that this routine is unable to recognise, then this routine will return
    the empty set.
    """
    if not arguments:
        print( "Need to supply at least one of the available test names:" )
        print( "    all" )
        for name in availableTests:
            print( "    {}".format(name) )
        print()
        return set()
    
    # At least one argument supplied. Check that all of these correspond to
    # known test names.
    testNames = set(arguments)
    unknown = testNames - set(availableTests) - {"all"}
    if unknown:
        print( "Unknown test name(s):" )
        for u in unknown:
            print( "    {}".format(u) )
        print( "Here are all the available test names:" )
        print( "    all" )
        for name in availableTests:
            print( "    {}".format(name) )
        print()
        return set()

    # All arguments are valid, so gather them all into a set of test names.
    if "all" in testNames:
        testNames = set(availableTests)
    print( "Running the following tests:" )
    for name in testNames:
        print( "    {}".format(name) )
    print()
    return testNames


def doTest( description, expected, actual ):
    """
    Tests that an actual computed value is equal to the expected value, and
    raises an AssertionError if this fails.
    """
    print( "{} Expected: {}. Actual: {}.".format(
        description, expected, actual ) )
    if expected != actual:
        raise AssertionError("TEST FAILED!")
    return


def allTestsPassedMessage(testNames):
    """
    Prints a message indicating that all the given tests passed.
    """
    if testNames:
        print( "+-------------------+" )
        print( "| ALL TESTS PASSED! |" )
        print( "+-------------------+" )
        for name in testNames:
            print(name)
