'''This module implements a class for retrieving, processing, and analyzing
signals from the interferometer channel of the mitpci system.

'''


# Standard library imports
import numpy as np
import matplotlib.pyplot as plt

# Related 3rd-party imports
from ..signal import Signal
from .ellipse import FittedEllipse


class Lissajous(object):
    '''An object containing the interferometer I&Q signals;
    plotting Q vs. I results in a Lissajous figure.

    Attributes:
    -----------
    shot - int
        The shot number of the retrieved signals.

    I - :py:class:`Signal <mitpci.signal.Signal>` instance
        The in-phase (I) signal of the interferometer.

    Q - :py:class:`Signal <mitpci.signal.Signal>` instance
        The quadrature (Q) signal of the interferometer.

    fit - bool
        If True, the I&Q signals have been fit to an ellipse
        (in the least-squares sense) and `self.E` will
        *not* be `None`.

    compensate - bool
        If True, the I&Q signals contained in `self.I.x` and
        `self.Q.x`, respectively, have been mapped from an
        ellipse to a circle using the fitting parameters
        in `self.E`.

    E - :py:class:`FittedEllipse
            <mitpci.interferometer.ellipse.FittedEllipse>` instance
        The object containing the elliptical parameters that
        result from least-squares fitting the raw I&Q signals.

        It is assumed that the properties of the ellipse
        vary adiabatically throughout the shot; that is,
        the properties of the ellipse vary much more slowly
        than the time taken for the raw I&Q signals
        to complete one full revolution around the origin.
        Each full revolution around the origin is then
        individually fit to an ellipse; compensation is
        similarly handled on a revolution-by-revolution basis.
        Because it is a time-consuming, iterative problem
        to both fit each full revolution *and* ensure
        continuity at the boundaries of each full revolution,
        continuity is *not* guaranteed at the revolution
        boundaries in the compensated I&Q signals. Instead,
        an approximation is used to guarantee continuity
        in the computed phase at the revolution boundaries,
        as the phase is what we ultimately care about.

    '''
    def __init__(
            self, shot, channel_I=1, channel_Q=2,
            fit=True, compensate=True,
            quiet=False, **signal_kwargs):
        '''Create an instance of the `Lissajous` class.

        Input parameters:
        -----------------
        shot - int
            The DIII-D shot number.

        channel_I - int
            The channel of the mitpci system corresponding to the
            interferometer's in-phase (I) signal

        channel_Q - int
            The channel of the mitpci system corresponding to the
            interferometer's quadrature (Q) signal

        fit - bool
            If True, fit the I&Q signals to an ellipse.
            If False but `compensate` is True, the I&Q signals
            are still fit, as the compensation depends on
            the fitting parameters.

        compensate - bool
            If True, use the fitted elliptical parameters to
            compensate the I&Q signals, mapping them from
            an ellipse to a circle.

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

        # Load raw I&Q signals
        if not quiet:
            print '\nRetrieving in-phase (I) signal for %i' % shot
        self.I = Signal(shot, channel_I, **signal_kwargs)
        self.I.x = self.I.x * self.I.volts_per_bit

        if not quiet:
            print 'Retrieving quadrature (Q) signal for %i' % shot
        self.Q = Signal(shot, channel_Q, **signal_kwargs)
        self.Q.x = self.Q.x * self.Q.volts_per_bit

        # Parse whether or not to fit and compensate I&Q
        if compensate:
            self.fit = True
            self.compensate = True
        elif fit:
            self.fit = True
            self.compensate = False
        else:
            self.fit = False
            self.compensate = False

        # Fit I&Q, if desired
        if self.fit:
            print '\nFitting I&Q signals to an ellipse'

            # `self.getFringeIndices()` returns indices corresponding
            # to the first point following a 2 * pi evolution
            # of the measured bulk phase; however, for reasons of
            # semantics, it does not return an index corresponding
            # to the initial data point. Manually include index
            # of initial data point here.
            starts = np.concatenate(([0], self.getFringeIndices()))

            self.E = FittedEllipse(
                self.I.x,
                self.Q.x,
                starts,
                t=self.I.t())
        else:
            self.E = None

        # Compensate I&Q, if desired
        if self.compensate:
            print '\nCompensating I&Q signals'
            self.I.x, self.Q.x = self.E.compensateEllipticity(
                self.I.x, self.Q.x)

    def getPhase(self, unwrap=True):
        'Return (potentially) unwrapped phase computed from I&Q signals.'
        # In our implementation, the LO power > the RF power.
        # For our I&Q demodulator (Mini-Circuits MIQC-60WD+),
        # LO > RF implies that the signal from the Q-pin obeys
        #
        #               Q(phi) = I(phi - pi / 2)
        #
        # That is, the signal from the Q pin lags the signal
        # from the I pin by pi / 2 radians. Now, if we assume
        # that I = I0 * cos(phi), this implies that
        # Q = Q0 * cos(phi - pi / 2) = Q0 * sin(phi), and
        #
        #               phi = tan^{-1}(Q / I)
        #
        # Note that `np.arctan2(...)` operates correctly
        # on integers and arrays of integers, so there is
        # no need to precondition `Q` and `I` as arrays
        # of floats
        ph = np.arctan2(self.Q.x, self.I.x)

        if unwrap:
            ph = np.unwrap(ph)

        return ph

    def getFringeIndices(self):
        '''Get indices corresponding to movement from one fringe to another
        (i.e. when the measured bulk phase evolves by 2 * pi).

        '''
        ph = self.getPhase()

        fringe_indices = _secular_change_indices(
            ph, secular_change=(2 * np.pi))

        return fringe_indices


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
            compensation={'DC': True, 'amplitude': True},
            quiet=False, **signal_kwargs):
        '''Create an instance of the `Demodulated` class.

        Input parameters:
        -----------------
        shot - int
            The shot number of the signal to be retrieved.

        channel_I [channel Q] - int
            The channel of the mitpci system corresponding to the
            interferometer's in-phase (I) [quadrature (Q)] signal

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

        self.compensation = compensation

        if self.compensation['DC']:
            print '\nSubtracting DC offsets from I&Q signals'
            self.subtractDCOffsets()

        if self.compensation['amplitude']:
            print 'Normalizing amplitude of Q to that of I'
            self.normalizeAmplitudes()

        return

    def getPhase(self, unwrap=True):
        'Return (potentially) unwrapped phase computed from I&Q signals.'
        # In our implementation, the LO power > the RF power.
        # For our I&Q demodulator (Mini-Circuits MIQC-60WD+),
        # LO > RF implies that the signal from the Q-pin obeys
        #
        #               Q(phi) = I(phi - pi / 2)
        #
        # That is, the signal from the Q pin lags the signal
        # from the I pin by pi / 2 radians. Now, if we assume
        # that I = I0 * cos(phi), this implies that
        # Q = Q0 * cos(phi - pi / 2) = Q0 * sin(phi), and
        #
        #               phi = tan^{-1}(Q / I)
        #
        # Note that `np.arctan2(...)` operates correctly
        # on integers and arrays of integers, so there is
        # no need to precondition `Q` and `I` as arrays
        # of floats
        ph = np.arctan2(self.Q.x, self.I.x)

        if unwrap:
            ph = np.unwrap(ph)

        return ph

    def getFringeIndices(self):
        '''Determine when interferometer passes from one fringe to another
        (i.e. when the measured phase evolves by 2 * pi).

        The indices are returned by the function call; additionally,
        the function assigns the indices to class attribute `fringe_indices`.

        '''
        ph = self.getPhase()

        self.fringe_indices = _secular_change_indices(
            ph, secular_change=(2 * np.pi))

        return self.fringe_indices

    def _getFringeSlice(self, fringe, fringe_indices):
        'Get slice corresponding to `fringe`.'
        if fringe == 0:
            sl = slice(0, fringe_indices[0])
        elif fringe == (len(fringe_indices) - 1):
            sl = slice(fringe_indices[-1], None)
        else:
            sl = slice(fringe_indices[fringe - 1], fringe_indices[fringe])

        return sl

    def _getFringeMax(self, x):
        'Get maximum value of `x` within each fringe.'
        try:
            fringe_indices = self.fringe_indices
        except AttributeError:
            fringe_indices = self.getFringeIndices()

        N = len(fringe_indices)
        fringe_max = np.zeros(N)

        for i in np.arange(N):
            sl = self._getFringeSlice(i, fringe_indices)
            fringe_max[i] = np.max(x[sl])

        return fringe_max

    def getDCOffsets(self):
        'Get I&Q DC offsets within each fringe.'
        # The I&Q signals may have been manipulated following any previous
        # computation of `self.fringe_indices`, so it is a good idea
        # to recompute the indices before performing further computations
        self.getFringeIndices()

        Imax = self._getFringeMax(self.I.x)
        Qmax = self._getFringeMax(self.Q.x)

        Imin = self._getFringeMax(-self.I.x)
        Qmin = self._getFringeMax(-self.Q.x)
        Imin *= -1
        Qmin *= -1

        # Ideally, the I&Q signals are pure sinusoids. If, however,
        # they have finite DC offset `dA`, their functional form becomes
        #
        #                   y = [A * cos(ph)] + dA
        #
        # As we are analyzing full fringes (i.e. `ph` passes through
        # a full 2 * pi radian), we know that
        #
        #                       ymax = A + dA
        #                       ymin = -A + dA
        #
        # and
        #
        #                   dA = (ymax + ymin) / 2
        dtype = self.I.x.dtype.name
        I_DC = (0.5 * (Imax + Imin)).astype(dtype)
        Q_DC = (0.5 * (Qmax + Qmin)).astype(dtype)

        return I_DC, Q_DC

    def subtractDCOffsets(self):
        'Subtract DC offsets within each fringe from I&Q signals.'
        I_DC, Q_DC = self.getDCOffsets()

        fringe_indices = self.fringe_indices
        N = len(fringe_indices)

        for i in np.arange(N):
            sl = self._getFringeSlice(i, fringe_indices)
            self.I.x[sl] -= I_DC[i]
            self.Q.x[sl] -= Q_DC[i]

        # As the I&Q signals have just been manipulated, it is a good idea
        # to recompute the indices before performing further computations
        self.getFringeIndices()

        self.compensation['DC'] = True

        return

    def getAmplitudes(self):
        'Get I&Q signal amplitudes within each fringe.'
        # The I&Q signals may have been manipulated following any previous
        # computation of `self.fringe_indices`, so it is a good idea
        # to recompute the indices before performing further computations
        self.getFringeIndices()

        Imax = self._getFringeMax(self.I.x)
        Qmax = self._getFringeMax(self.Q.x)

        Imin = self._getFringeMax(-self.I.x)
        Qmin = self._getFringeMax(-self.Q.x)
        Imin *= -1
        Qmin *= -1

        # Assuming a pure sinusoid with a DC offset, the I&Q signals
        # have the form
        #
        #                   y = [A * cos(ph)] + dA
        #
        # As we are analyzing full fringes (i.e. `ph` passes through
        # a full 2 * pi radian), we know that
        #
        #                       ymax = A + dA
        #                       ymin = -A + dA
        #
        # and
        #
        #                   A = (ymax - ymin) / 2
        dtype = self.I.x.dtype.name
        I0 = (0.5 * (Imax - Imin)).astype(dtype)
        Q0 = (0.5 * (Qmax - Qmin)).astype(dtype)

        return I0, Q0

    def normalizeAmplitudes(self):
        'Normalize amplitude of Q to that of I within each fringe.'
        I0, Q0 = self.getAmplitudes()

        Q = self.Q
        dtype = self.Q.x.dtype.name

        fringe_indices = self.fringe_indices
        N = len(fringe_indices)

        for i in np.arange(N):
            sl = self._getFringeSlice(i, fringe_indices)
            norm = np.float(I0[i]) / Q0[i]
            Q.x[sl] = (Q.x[sl] * norm).astype(dtype)

        # As the I&Q signals have just been manipulated, it is a good idea
        # to recompute the indices before performing further computations
        self.getFringeIndices()

        self.compensation['amplitude'] = True

        return


def _secular_change_indices(x, secular_change=1, plot=False):
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
