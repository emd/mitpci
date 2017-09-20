'''This module implements a class for retrieving and processing
PCI signals digitized on the mitpci system.

'''


# Standard library imports
import numpy as np
import matplotlib.pyplot as plt

# Related 3rd-party imports
import filters
from random_data.signals.sampling import circular_resample

# Intra-package imports
from ..signal import Signal
from ..interferometer.demodulated import _hpf
from .detector_array import DetectorArray


class Phase(Signal):
    '''An object containing the PCI-measured phase, potentially
    from several channels. Note that the phase is only determined
    up to a calibration constant.

    This class is derived from :py:class:`Signal <mitpci.signal.Signal>`
    and thus shares most of its attributes and methods. Attributes and
    properties *unique* to :py:class:`Phase <mitpci.pci.signal.Phase>`
    are discussed below.

    Attributes:
    -----------
    shot - int
        The shot number of the phase signal.

    x - array-like, (`M`, `N`) or, if only one channel, (`N`,)
        The PCI-measured phase, determined up to a calibration constant.
        Here `M` is the number of PCI channels and `N` is the number
        of timestamps in the retrieved signal.
        [x] = rad

    Fs - float
        The signal sampling rate.
        [Fs] = samples / second

    t0 - float
        The time corresponding to `self.x[0]`.
        [t0] = s

    digitizer_channels - array_like, (`M`,)
        The digitizer channels corresponding to `self.x`. That is,
        `self.x[i, :]` corresponds to the signal retrieved from
        digitizer channel `self.digitizer_channels[i]`.
        [digitizer_channels] = unitless

    digitizer_board - array_like, (`M`,)
        The digitizer board for each channel in `digitizer_channels`.
        [digitizer_board] = unitless

    detector_elements - array_like, (`M`,)
        The detector elements corresponding to `self.x`. That is,
        `self.x[i, :]` corresponds to the signal measured from
        digitizer element `self.detector_elements[i]`.
        [detector_elements] = unitless

    detector_stencil - :py:class:`DetectorArray
            <mitpci.pci.detector_array.DetectorArray`
        The stencil corresponding to `self.detector_elements`.

    filt - :py:class:`Kaiser <filters.fir.Kaiser>` instance or None
        The filter applied to the phase signal. If `None`, the phase
        signal has not been filtered. Signal and time points
        contaminated by the filter's boundary effects are *not*
        returned/accessible from the attributes and methods
        of the `Phase` class.

    '''
    def __init__(self, shot, digitizer_channels, tau=None,
                 filt=_hpf, quiet=False, **signal_kwargs):
        '''Create an instance of the `Phase` class.

        Input parameters:
        -----------------
        shot - int
            The DIII-D shot number.

        digitizer_channels - array_like, (`L`,)
            The digitizer channels from which to retrieve PCI data.
            This may be the full set of available digitizer channels
            or any subset thereof. Digitizer channels corresponding
            to dead detector elements or mapped to something other
            than the PCI detector array (e.g. the I&Q signals of the
            heterodyne interferometer) are automatically removed. A
            ValueError is raised if non-existent digitizer channels are
            specified.

        tau - float or None
            The trigger offset between boards 7 and 8. If the board-7
            measurements correspond to a true, physical time `t`, then
            the board-8 measurements correspond to a true, physical
            time `t + tau`. (Note that `tau` is *signed* such that a
            negative value of `tau` indicates that the true, physical
            time corresponding to the board-8 measurements actually
            precedes that of the board 7 measurements).

            If `tau` is not `None`, then the trigger offset will be
            compensated. If `tau` is `None` and data from both boards
            are loaded, a cautionary note is printed to the screen.

            The trigger offset can be easily estimated using
            py:class:`TriggerOffset <random_data.signals.TriggerOffset>`.

            [tau] = s

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

                    Ph = Phase(167340, [7, 8], tlim=[1, 3.3])

            to retrieve the phase signals measured by digitizer
            channels 7 & 8 of the PCI system for shot 167340 between
            1 <= t [s] <= 3.3.

        '''
        # If single channel passed as an integer, convert to `array_like`
        if type(digitizer_channels) is int:
            digitizer_channels = np.array([digitizer_channels])

        self.tau = tau

        # Create stencil from specified digitizer channels
        self.detector_stencil = DetectorArray(shot, digitizer_channels)

        # Load raw signal, compensating for the trigger offset
        # between board 7 and board 8 if needed
        self._getRawSignal(quiet=quiet, **signal_kwargs)

        # Filter signal, if requested
        self.filt = filt
        if self.filt is not None:
            self._filterSignal()

        # Convert from bits to radians using measurements
        # from PCI and interferometer cross calibration
        # with swept-frequency sound waves
        rad_per_bit = 1. / 7.7e5
        self.x = self.x * rad_per_bit

        Nboards = len(set(self.digitizer_board))

        # If trigger offset is *not* compensated but signals from
        # both boards have been loaded, warn the user
        if (self.tau is None) and (Nboards > 1):
            print '\nWARNING: Trigger offset between boards NOT compensated,'
            print 'which may result in significant systematic errors.'
            print 'Offset estimation & specification discussed in docs.'

    def _getRawSignal(self, quiet=False, **signal_kwargs):
        'Get raw signal from `self.digitizer_channels`.'
        # Load raw signal from first channel
        if not quiet:
            self._printLoadingMessage(0)

        shot = self.detector_stencil.shot
        ch = self.digitizer_channels[0]
        Signal.__init__(self, shot, ch, **signal_kwargs)

        # Create a list to store digitizer board of channel(s)
        self.digitizer_board = [self._digitizer_board]

        # Delete useless & confusing attributes for this derived class
        del self.channel, self._digitizer_board, self._node_name

        Nchannels = len(self.digitizer_channels)

        if Nchannels > 1:
            # Initialize an array to hold data from all of the
            # specified digitizer channels
            Ntimes = len(self.x)
            xtmp = self.x
            xdtype = self.x.dtype
            self.x = np.zeros((Nchannels, Ntimes), dtype=xdtype)
            self.x[0, :] = xtmp

            # Loop through remaining channels
            for channel_index in (np.arange(Nchannels - 1) + 1):
                if not quiet:
                    self._printLoadingMessage(channel_index)

                ch = self.digitizer_channels[channel_index]
                sig = Signal(shot, ch, **signal_kwargs)
                self.x[channel_index, :] = sig.x
                self.digitizer_board.append(sig._digitizer_board)

        # Convert to array for compatibility with other arrays,
        # such as `self.digitizer_channels` and `self.detector_elements`
        self.digitizer_board = np.array(self.digitizer_board)
        Nboards = len(set(self.digitizer_board))

        # If needed, compensate for trigger offset
        if (self.tau is not None) and (Nboards > 1):
            print ''
            for chind in np.where(self.digitizer_board == 'DT216_8')[0]:
                self._printOffsetCompensationMessage(chind)
                self.x[chind, :] = circular_resample(
                    self.x[chind, :],
                    self.Fs,
                    -self.tau,
                    pad_to_next_fast_len=True)

        return

    def _printLoadingMessage(self, channel_index):
        'Print digitizer channel number and shot when loading data.'
        ch = self.detector_stencil.digitizer_channels[channel_index]
        shot = self.detector_stencil.shot

        if channel_index == 0:
            prepend = '\n'
        else:
            prepend = ''

        print '%sRetrieving PCI ch. %i for %i' % (prepend, ch, shot)

        return

    def _printOffsetCompensationMessage(self, channel_index):
        'Print channel number when compensating for trigger offset.'
        ch = self.digitizer_channels[channel_index]
        tau = self.tau * 1e6

        print 'Compensating PCI ch. %i for %.2f us offset' % (ch, tau)

        return

    def _filterSignal(self):
        'Filter signal `self.x` with `self.filt`.'
        # Ensure `filt` is of correct type
        if isinstance(self.filt, filters.fir.Kaiser):
            if self.filt.Fs == self.Fs:
                valid = self.filt.getValidSlice()
            else:
                raise ValueError(
                    '`filt` not designed for signal sample rate')
        else:
            raise ValueError(
                '`filt` must be of type `filters.fir.Kaiser`')

        Nchannels = len(self.digitizer_channels)

        if Nchannels == 1:
            # Apply filter to single channel. As `self.x` is only
            # 1-dimensional in the case of a single channel, we
            # must use different syntax in the filter application.
            print ('\nFiltering PCI ch. %i for %i' %
                   (self.digitizer_channels[0], self.shot))
            self.x = self.filt.applyTo(self.x)
        else:
            # Loop through each channel to apply filter
            for ind in np.arange(Nchannels):
                print ('\nFiltering PCI ch. %i for %i' %
                       (self.digitizer_channels[ind], self.shot))
                self.x[ind, :] = self.filt.applyTo(self.x[ind, :])

        # Only allow points that are *free* from the filter's boundary
        # effects to be accessible; modify initial time accordingly
        self.x = self.x[..., valid]
        self.t0 = self.t()[valid][0]

        return

    @property
    def digitizer_channels(self):
        'Get digitizer channels corresponding to `self.x`.'
        return self.detector_stencil.digitizer_channels

    @property
    def detector_elements(self):
        'Get detector elements corresponding to `self.x`.'
        return self.detector_stencil.detector_elements
