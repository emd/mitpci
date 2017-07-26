# Standard library imports
import numpy as np
import matplotlib.pyplot as plt

# Related 3rd-party imports
import filters

# Intra-package imports
from ..signal import Signal
from ..interferometer.demodulated import _hpf


class Phase(Signal):
    '''An object containing the PCI-measured phase. Note that
    the phase is only determined up to a calibration constant.

    Attributes:
    -----------
    shot - int
        The shot number of the phase signal.

    channel - int
        The channel of the PCI system.

    filt - :py:class:`Kaiser <filters.fir.Kaiser>` instance or None
        The filter applied to the phase signal. If `None`, the phase
        signal has not been filtered. Signal and time points
        contaminated by the filter's boundary effects are *not*
        returned/accessible from the attributes and methods
        of the `Phase` class.

    x - array-like, (`N`,)
        The PCI-measured phase, determined up to a calibration constant.
        [x] = rad

    Fs - float
        The signal sampling rate.
        [Fs] = samples / second

    t0 - float
        The time corresponding to `self.x[0]`.
        [t0] = s

    '''
    def __init__(self, shot, channel, filt=_hpf,
                 quiet=False, **signal_kwargs):
        '''Create an instance of the `Phase` class.

        Input parameters:
        -----------------
        shot - int
            The DIII-D shot number.

        channel - int
            The channel of the PCI system.

        filt - :py:class:`Kaiser <filters.fir.Kaiser>` instance or None
            The filter applied to the phase signal. If `None`,
            do not filter the phase signal. A ValueError is raised
            if `filt` is *not* an instance of `filters.fir.Kaiser`
            or `None`.

        quiet - bool
            If True, suppress printing messages to the terminal.

        signal_kwargs - any valid keyword arguments for
            :py:class:`Signal <mitpci.signal.Signal>`.

            For example, use

                    Ph = Phase(167340, 8, tlim=[1, 3.3])

            to retrieve the phase signal measured by channel 8
            of the PCI system for shot 167340 between
            1 <= t [s] <= 3.3.

        '''
        if not quiet:
            print '\nRetrieving PCI ch. %i for %i' % (channel, shot)

        # Load raw signal
        Signal.__init__(
            self, shot, channel, **signal_kwargs)

        # Ensure `filt` is of correct type
        if isinstance(filt, filters.fir.Kaiser):
            if filt.Fs == self.Fs:
                self.filt = filt
                self._valid = self.filt.getValidSlice()
            else:
                raise ValueError('`filt` not designed for signal sample rate')
        elif filt is None:
            self.filt = None
            self._valid = slice(None, None)
        else:
            raise ValueError(
                '`filt` must be `filters.fir.Kaiser` or `None`')

        # Apply filter, returning only points that are free from
        # the filter's boundary effects; modify initial time accordingly
        if self.filt is not None:
            self.x = self.filt.applyTo(self.x)[self._valid]
            self.t0 = self.t()[self._valid][0]

        # Convert from bits to radians using measurements
        # from PCI and interferometer cross calibration
        # with swept-frequency sound waves
        rad_per_bit = 1. / 7.7e5
        self.x = self.x * rad_per_bit
