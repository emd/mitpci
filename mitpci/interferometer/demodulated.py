'''This module implements a class for retrieving, processing, and analyzing
signals from the interferometer channel of the mitpci system.

'''


# Standard library imports
import numpy as np
import matplotlib.pyplot as plt

# Related 3rd-party imports
from ..signal import Signal
from .ellipse import FittedEllipse
import filters


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
        If `fit` is False, `E` will be `None`.

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
        self.fit = False
        self.compensate = False
        self.E = None

        # Load raw I&Q signals
        if not quiet:
            print '\nRetrieving in-phase (I) signal for %i' % shot
        self.I = Signal(shot, channel_I, **signal_kwargs)
        self.I.x = self.I.x * self.I.volts_per_bit

        if not quiet:
            print 'Retrieving quadrature (Q) signal for %i' % shot
        self.Q = Signal(shot, channel_Q, **signal_kwargs)
        self.Q.x = self.Q.x * self.Q.volts_per_bit

        # Fit and compensate I&Q, if desired
        if fit or compensate:
            self.getFit(quiet=quiet)

            if compensate:
                self.compensateEllipticity(quiet=quiet)

    def getPhase(self, unwrap=True, enforce_fringe_continuity=True):
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

        # Compensation of I&Q ellipticity on a fringe-by-fringe basis
        # does *not* ensure continuity of the resulting phase
        # across fringe boundaries. Here, if desired, we force
        # the phase change across each fringe boundary to be equal
        # to the corresponding phase change computed from the
        # *raw* I&Q signals, ensuring continuity across fringes.
        # Note that this may result in an offset in the bulk phase;
        # however, the bulk phase in our system corresponds to
        # vibrations and is of little physical importance --
        # we only care that the small, plasma-induced fluctuations
        # that ride on top of the vibration-induced phase
        # are mapped onto a circle.
        if self.compensate and enforce_fringe_continuity:
            ph = _enforce_boundary_conditions(
                ph, self.E.starts, self._raw_phase_change_at_fringe)

        return ph

    def getFit(self, quiet=False):
        '''Get elliptical fit for each revolution of `self.I.x` and
        `self.Q.x` about origin.

        '''
        # I&Q signals *not* already fit
        if not self.fit:
            if not quiet:
                print '\nFitting I&Q signals to an ellipse'

            fringe_info = self._getRawFringeInfo(quiet=quiet)
            starts = fringe_info[0]
            self._raw_phase_change_at_fringe = fringe_info[1]

            self.E = FittedEllipse(
                self.I.x,
                self.Q.x,
                starts,
                t=self.I.t())

            self.fit = True

        # I&Q signals already fit
        else:
            if not quiet:
                print '\nI&Q signals already fit to an ellipse'

        # `self.E` is not really intended for use outside
        # of object instance, so don't return anything
        return

    def compensateEllipticity(self, quiet=False):
        '''Use fitting parameters to compensate ellipticity
        of measured I&Q signals, effectively mapping
        `self.I.x` and `self.Q.x` onto a circle.

        '''
        # I&Q signals *not* already compensated
        if not self.compensate:
            # Fit ellipse, if not already fit
            if not self.fit:
                self.getFit(quiet=quiet)

            if not quiet:
                print '\nCompensating I&Q signals'

            self.I.x, self.Q.x = self.E.compensateEllipticity(
                self.I.x, self.Q.x)

            self.compensate = True

        # I&Q signals already compensated
        else:
            if not quiet:
                print '\nI&Q signals already compensated'

        return

    def _getRawFringeInfo(self, quiet=False):
        '''Get indices corresponding to the start of a new fringe
        in the *raw* I&Q signals (i.e. when the raw, measured bulk
        phase evolves by 2 * pi) and the change in phase at the
        fringe boundary.

        '''
        # I&Q signals *not* already fit (i.e. need to compute fringe info)
        if not self.fit:
            ph = self.getPhase()

            ind = _secular_change_indices(ph, secular_change=(2 * np.pi))

            # `_secular_change_indices` returns indices corresponding
            # the first point following each full revolution of the I&Q
            # signals about the origin; however, for reasons of
            # semantics, it does not return an index corresponding
            # to the initial data point. Manually include index
            # of initial data point here.
            ind = np.concatenate(([0], ind))

            # Determine change in phase at boundary between
            # successive fringes. Note that there is *not*
            # a fringe prior to `ind[0]`; thus, `ind[0]` does
            # *not* correspond to a boundary between two fringes.
            dph = _get_boundary_delta(ph, ind)

            return ind, dph

        # I&Q signals already fit (i.e. do not need to compute fringe info)
        else:
            if not quiet:
                print '\nRaw fringe information already calculated'

            return


class Phase(object):
    '''An object containing the interferometer-measured phase.

    Attributes:
    -----------
    shot - int
        The shot number of the phase signal.

    compensate - bool
        If True, systematic errors in the phase calculation
        have been compensated/minimized by mapping the I&Q signals
        (from which the phase is derived) from an ellipse to a circle.

    filt - :py:class:`Kaiser <filters.fir.Kaiser>` instance or None
        The filter applied to the phase signal. If `None`, the phase
        signal has not been filtered. Signal and time points
        contaminated by the filter's boundary effects are *not*
        returned/accessible from the attributes and methods
        of the `Phase` class.

    x - array-like, (`N`,)
        The retrieved phase signal, as derived from the I&Q signals
        of input `L` and filtered with `self.filt`. Only points free
        of the filter's boundary effects are accessible.
        [x] = rad

    Fs - float
        The signal sampling rate.
        [Fs] = samples / second

    t0 - float
        The time corresponding to `self.x[0]`; note that `self.x[0]`
        and `self.x[-1]` (and all the points in between) are *free*
        from the boundary effects resulting from application of
        `self.filt`.
        [t0] = s

    '''
    def __init__(self, L, filt=None):
        '''Create an instance of the `Phase` class.

        Parameters:
        -----------
        L - :py:class:`Lissajous
                <mitpci.interferometer.demodulated.Lissajous>` instance
            Instance of `Lissajous` class containing I&Q signals
            for given shot. A ValueError is raised if `L` is *not*
            an instance of the `Lissajous` class.

        filt - :py:class:`Kaiser <filters.fir.Kaiser>` instance or None
            The filter applied to the phase signal. If `None`,
            do not filter the phase signal. A ValueError is raised
            if `filt` is *not* an instance of `filters.fir.Kaiser`
            or `None`.

        '''
        # Ensure `L` is of correct type
        if not isinstance(L, Lissajous):
            raise ValueError(
                '`L` must be `mitpci.interferometer.demodulated.Lissajous`')

        self.shot = L.shot
        self.compensate = L.compensate
        self.Fs = L.I.Fs

        # Ensure `filt` is of correct type
        if isinstance(filt, filters.fir.Kaiser):
            self.filt = filt
            self._valid = self.filt.getValidSlice()
        elif filt is None:
            self.filt = None
            self._valid = slice(None, None)
        else:
            raise ValueError(
                '`filt` must be `filters.fir.Kaiser` or `None`')

        self.x = self._getPhase(L)
        self.t0 = L.I.t()[self._valid][0]

    def _getPhase(self, L):
        'Get phase signal after filtering by `self.filt`.'
        ph = L.getPhase()

        if self.filt is not None:
            ph = self.filt.applyTo(ph)

        return ph[self._valid]

    def t(self):
        'Get times for points in `self.x`.'
        return self.t0 + (np.arange(len(self.x)) / self.Fs)


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


def _get_boundary_delta(x, starts):
    '''Get change in `x` between `starts` and `starts - 1`.

    Parameters:
    -----------
    x - array_like, (`N`,)
        An array that has been subdivided into `M` sub-segments.

    starts - array_like, (`M`,)
        Index specifying the start of each sub-segment in `x`,
        with the ith sub-segment beginning at index `starts[i]`.
        `M` < `N`

    Returns:
    --------
    dx - array_like, (`M` - 1,)
        Array corresponding to the *change* in `x` across
        the boundary of each sub-segment.

    '''
    if len(starts) > 1:
        dx = x[starts[1:]] - x[starts[1:] - 1]
    else:
        # Only one sub-segment, so there is, by definition,
        # no boundary between successive sub-segments
        dx = None

    return dx


def _enforce_boundary_conditions(x, starts, bc):
    '''Get array corresponding to `x` altered by boundary conditions `bc`.

    Parameters:
    -----------
    x - array_like, (`N`,)
        An array that has been subdivided into `M` sub-segments.
        The change in `x` between the ith and (i - 1)th segment
        will be changed to the boundary condition specified
        by `bc[i]`.

    starts - array_like, (`M`,)
        Index specifying the start of each sub-segment in `x`,
        with the ith sub-segment beginning at index `starts[i]`.
        `M` < `N`

    bc - array_like, (`M` - 1,)
        The desired boundary conditions between successive sub-segments.
        Denoting the returned, corrected signal as `xc`, the boundary
        condition between the ith and (i - 1)th segment of `xc` is

            bc[i] = xc[starts[i]] - xc[starts[i] - 1]

    Returns:
    --------
    xc - array_like, (`N`,)
        The array corresponding to input signal `x` modified by
        the boundary conditions specified in `bc`.

    '''
    xc = x.copy()

    # Get existing boundary conditions of input array `x`
    dx = _get_boundary_delta(x, starts)

    # To enforce the desired boundary conditions such that
    #
    #       bc[i] = xc[starts[i]] - xc[starts[i] - 1],
    #
    # we must *subtract* the existing boundary condition `dx[i]`
    # from all points beyond the (i - 1)th segment and
    # *add* the desired boundary conditions `bc[i]`
    # to all points beyond the (i - 1)th segment.
    #
    # Perhaps the most straightforward approach is iterative; i.e.
    #
    #       for i, start in enumerate(starts):
    #           xc[start:] -= (dx[i] - bc[i])
    #
    # However, for large `x` (e.g. millions of data points) and
    # large `starts` (e.g. splitting `x` into thousands of segments),
    # the above iterative approach is *very slow*.
    #
    # Note that the above boundary-condition enforcement
    # for the ith boundary is equivalent to a *cumulative sum*
    # of all of the preceding boundary conditions up to and
    # including the ith boundary. The computation can thus
    # be partially vectorized and greatly accelerated
    # by using Numpy's `cumsum(...)` function.
    dx_cum = np.cumsum(dx)
    bc_cum = np.cumsum(bc)

    Nbnd = len(starts[1:])

    for bnd in np.arange(Nbnd):
        # Initial index for slicing segment
        i1 = starts[1:][bnd]

        # Final index for slicing segment, with
        # appropriate handling of last segment
        if bnd < (Nbnd - 1):
            i2 = starts[1:][bnd + 1]
        else:
            i2 = None

        # Slice segment and enforce boundary condition
        sl = slice(i1, i2)
        xc[sl] -= (dx_cum[bnd] - bc_cum[bnd])

    return xc
