import numpy as np
from mitpci.interferometer.demodulated import _secular_change_indices


def test_secular_change_indices():
    # With the below definition for the sawtooth function, we see that
    # for points `x1` and `x2` satisfying floor(x1) == floor(x2)
    #
    #               x2 - x1 = (y2 - y1) / (2 * A)
    #
    # Now, to test our `_secular_change_indices(...)` function,
    # which has a default secular change of unity (i.e. y2 - y1 = 1),
    # we should choose the spacing of our computational grid to be
    #
    #                       dx = 1 / (2 * A)
    #
    # such that each grid point is identified as having a "secular change".
    # Of course, as discussed in the source for `_secular_change_indices(...)`,
    # the first "jump" does not have a preceding "jump" to be compared against,
    # so we expect grid points 0 and 1 to *not* have a secular change.

    A = 4  # to avoid round-off errors, A *needs* to be a power of 2!
    dx = 1. / (2 * A)

    num_cycles = 5
    x = np.arange(0, num_cycles, dx)

    y = A * sawtooth(x)
    ind = _secular_change_indices(y)

    ind_expected = np.arange(2, (2 * A * num_cycles))

    np.testing.assert_equal(ind, ind_expected)

    return


def sawtooth(x):
    '''Return a rising sawtooth with range [-1, 1] and unity period.
    (The Scipy sawtooth function is not working quite as expected, and
    it is making testing difficult; thus, we have this simpler function).

    '''
    return (2 * (x - np.floor(x))) - 1
