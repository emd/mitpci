import numpy as np
from nose import tools
import filters
from mitpci.interferometer.demodulated import (
    _secular_change_indices, _get_boundary_delta,
    _enforce_boundary_conditions,
    Lissajous, Phase)


def test__secular_change_indices():
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


def test__get_boundary_delta():
    da = 10
    a = np.arange(0, 110, da)
    starts = np.array([0, 5])

    np.testing.assert_equal(
        np.array([da]),
        _get_boundary_delta(a, starts))

    prefactor = 10
    np.testing.assert_equal(
        np.array([prefactor * da]),
        _get_boundary_delta(prefactor * a, starts))

    divisor = 10
    np.testing.assert_equal(
        np.array([da / divisor]),
        _get_boundary_delta(a / divisor, starts))

    return


def test__enforce_boundary_conditions():
    # Create initial array
    astart = 0
    astop = 110
    da = 1
    a = np.arange(astart, astop, da)

    # Specify sub-segments of `a` and
    # corresponding boundary conditions
    starts = np.arange(astart, astop, 10 * da)
    bc = da * np.ones(len(starts) - 1).astype(a.dtype.name)

    # Make copies of `a` whose boundary conditions
    # will subsequently be perturbed
    a1 = a.copy()
    a2 = a.copy()

    # Perturb boundary conditions of `a1` and `a2`
    for i, start in enumerate(starts[1:]):
        a1[start:] += (((-1) ** i) * start)
        a2[start:] += (((-1) ** (i + 1)) * start)

    # Enforcing original boundary conditions on `a1` and `a2`
    # should result in equality with `a`
    np.testing.assert_equal(
        a,
        _enforce_boundary_conditions(a1, starts, bc))

    np.testing.assert_equal(
        a,
        _enforce_boundary_conditions(a2, starts, bc))

    # Further, original boundary conditions on `a`
    # should *not* alter `a`
    np.testing.assert_equal(
        a,
        _enforce_boundary_conditions(a, starts, bc))

    return


def test_Phase__init__():
    shot = 171110
    L = Lissajous(shot, fit=False, compensate=False)

    # Needs to receive object of type `Lissajous`
    tools.assert_raises(
        ValueError,
        Phase, *[L.I], **{'filt': None})

    # Needs to receive `filt` that is `None` or `filters.fir.Kaiser`
    tools.assert_raises(
        ValueError,
        Phase, *[L], **{'filt': {}})

    # Use incorrect sampling rate
    Fs_bad = 0.5 * L.I.Fs
    hpf = filters.fir.Kaiser(-120, 5e3, 10e3, pass_zero=False, Fs=Fs_bad)
    tools.assert_raises(
        ValueError,
        Phase, *[L], **{'filt': hpf})

    return
