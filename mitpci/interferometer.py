'''This module implements a class for retrieving, processing, and analyzing
signals from the interferometer channel of the mitpci system.

'''


# Standard library imports
import numpy as np
import matplotlib.pyplot as plt

# Related 3rd-party imports
from .signal import Signal


class Demodulated(object):
    '''An object corresponding to the I&Q signals of the interferometer.

    Attributes:
    -----------
    shot - int
        The shot number of the retrieved signals.

    I - :py:class:`Signal <mitpci.signal.Signal>` instance
        The in-phase (I) signal of the interferometer.

    Q - :py:class:`Signal <mitpci.signal.Signal>` instance
        The quadrature (Q) signal of the interferometer.

    Methods:
    --------

    '''
    def __init__(
            self, shot, channel_I=1, channel_Q=2,
            quiet=False, **signal_kwargs):
        '''Create an instance of the `Demodulated` class.

        Input parameters:
        -----------------
        shot - int
            The shot number of the signal to be retrieved.

        channel_I [channel Q] - int
            The channel of the mitpci system corresponding to the
            interferometer's in-phase (I) [quadrature (Q) signal]signal

        quiet - bool
            If True, suppress printing messages to the terminal.

        signal_kwargs - any valid keyword arguments for
            :py:class:`Signal <mitpci.signal.Signal>`.

            For example, use

                    D = Demodulated(167340, tlim=[1, 3.3])

            to retrieve the I and Q signals for shot 167340 between
            1 <= t [s] <= 3.3 from channels 1 and 2, respectively,
            of the mitpci system.

        '''
        self.shot = shot

        if not quiet:
            print '\nRetrieving in-phase (I) signal for %i' % shot
        self.I = Signal(shot, channel_I, **signal_kwargs)

        if not quiet:
            print 'Retrieving quadrature (Q) signal for %i' % shot
        self.Q = Signal(shot, channel_Q, **signal_kwargs)

        return


def secular_change_indices(x, secular_change=1, plot=False):
    '''Get indices between which `x` changes by `secular_change`.

    Parameters:
    -----------
    x - array_like, (`N`,)
        The signal to be examined for secular variation.

    secular_change - float
        The minimum amount that `x` must change by to be considered
        a secular variation (as opposed to e.g. a small-amplitude
        periodic variation)

    plot - bool
        If True, plot signal `x` and its corresponding windows
        of secular change. Plotting is not recommended for long signals.

    Returns:
    --------
    ind - array_like, (`M`) with `M` < `N`
        The indices for which `x` changes by `secular_change`;
        that is, for `i` in {1, 2, ... `M`}, we will have:

            |x[ind[i]] - x[ind[i - 1]]| >= secular_change

        The boundary case for i = 0 satisfies

                |x[ind[0]] - x[0]| >= secular_change

    '''
    # Convert "analog" signal `x` into a "digital" signal
    # with minimum bit size `secular_change`
    xdig = np.floor(x / secular_change) * secular_change

    # Find positions where `xdig` jumps from one "bit" to another,
    # with `jump` storing the indices of the *first* point
    # following each jump
    jump = np.where(np.diff(xdig) != 0)[0] + 1

    # Handle the (rare) case that `jump[-1]` exceeds the bounds
    # of the signal
    if jump[-1] == len(xdig):
        jump = jump[:-1]

    # To be a *secular* jump (rather than a jump corresponding
    # to e.g. a small-amplitude periodic signal), the value of
    # the "digital" signal at `jump[i]` must *not* be equal to
    # its value immediately prior to `jump[i - 1]`. (Sketching
    # the  "analog" and "digital" signals may aid in the
    # understanding of this point).
    secular_jump = np.where(
        xdig[jump[1:]] != xdig[jump[:-1] - 1])[0]

    # The secular nature of any given jump is dependent
    # on the signal value immediately prior to the preceding jump.
    # Now, the *first* jump does *not* have a preceding jump, and
    # we can therefore not determine whether or not it is
    # a secular jump. This explains the `[1:]` indexing below.
    ind = jump[1:][secular_jump]

    if plot:
        plt.figure()

        plt.plot(x)
        plt.plot(xdig)
        plt.plot(jump, xdig[jump], 'ro')
        for i in ind:
            plt.axvline(i, c='k')

        plt.xlabel('index')
        plt.ylabel('value')
        plt.legend(
            ['signal', '"digital" signal', '"jumps"', 'secular changes'],
            loc='upper right')

        plt.show()

    return ind
