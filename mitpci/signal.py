'''This module implements a class for interacting with the signals
digitized by the `mitpci` computer system.

The `mitpci` computer system digitizes data from MIT's (a) phase contrast
imaging (PCI) and (b) heterodyne interferometer systems. The systems share
a beam path through the DIII-D vessel and a digitizer, but they use different
interference schemes and detectors to make their respective measurements.

'''


import numpy as np
import MDSplus as mds


class Signal(object):
    def __init__(self, shot, channel, channels_per_board=8,
                 Fs=None, tlim=None):
        self.shot = shot
        self.channel = channel

        # Check that `self.channel` is valid
        self._checkChannels(channels_per_board=channels_per_board)

        # Obtain digitizer board and signal node name for `self.channel`
        self._digitizer_board = self._getDigitizerBoard(
            channels_per_board=channels_per_board)
        self._node_name = self._getNodeName(
            channels_per_board=channels_per_board)

        # Open the tree and retrieve the signal within the specified
        # time window and with the specified sampling rate
        mds_tree = mds.Tree('pci', shot=shot, mode='ReadOnly')
        self.Fs, self._downsample = self._getSampleRate(mds_tree, Fs=Fs)
        if shot != -1:
            self.t0, self.x = self._getSignal(mds_tree, tlim=tlim)

    def _checkChannels(self, channels_per_board=8):
        'Ensure that `self.channel` corresponds to a physical mitpci channel.'
        if (self.channel <= 0) or (self.channel > (2 * channels_per_board)):
            raise ValueError(
                'Valid `channel` values between 1 and %i' %
                (2 * channels_per_board))
        return

    def _getDigitizerBoard(self, channels_per_board=8):
        'Get digitizer board corresponding to `self.channel`.'
        if self.channel <= channels_per_board:
            return 'DT216_7'
        else:
            return 'DT216_8'

    def _getNodeName(self, channels_per_board=8):
        board_channel = 1 + ((self.channel - 1) % channels_per_board)
        node_name = '.HARDWARE:%s:INPUT_%s' % (self._digitizer_board,
                                               str(board_channel).zfill(2))
        return node_name

    def _getSampleRate(self, mds_tree, Fs=None):
        'Get signal sampling rate, including effects of desired downsampling.'
        node = mds_tree.getNode(
            '.HARDWARE:%s:CLOCK_DIV' % self._digitizer_board)

        digitizer_rate = np.float(node.data())

        # Downsample if desired
        if Fs is not None:
            if Fs > 0:
                downsample = np.int(np.floor(digitizer_rate / Fs))
            else:
                raise ValueError('`Fs` must be positive')
        else:
            downsample = 1

        return (digitizer_rate / downsample, downsample)

    def _getSlice(self, x, tlim=None, t0_dig=0.):
        'Get slice for signal retrieval between `tlim` at rate `self.Fs`.'
        # Minimum and maximum values for slicing `x`
        imin = 0
        imax = len(x)

        if tlim is not None:
            # Ensure limits in time are correctly sized and sorted
            if len(tlim) != 2:
                raise ValueError('`tlim` must be an iterable of length 2.')
            else:
                tlim = np.sort(tlim)

            # Digitization rate
            Fs_dig = self.Fs * self._downsample

            # Find slicing indices such that:
            #   (a) `x[ilo:ihi]` corresponds to the signal within `tlim`, and
            #   (b) `ilo` and `ihi` are bounded by `imin` and `imax`
            #
            ilo = np.max([
                imin,
                np.ceil((tlim[0] - t0_dig) * Fs_dig)])

            # If `ihi_exact` is an integer, then `tlim[1]`
            # sits *exactly* on a digitized point; to ensure
            # we include this point in our slice, we should
            # augment `ihi_exact` by +1. If `ihi_exact` is
            # *not* an integer, the ceiling operation takes
            # care of this concern for us.
            ihi_exact = (tlim[1] - t0_dig) * Fs_dig
            if ihi_exact != np.int(ihi_exact):
                ihi = np.min([imax, np.ceil(ihi_exact)])
            else:
                ihi = np.min([imax, ihi_exact + 1])
        else:
            ilo = imin
            ihi = imax

        return slice(ilo, ihi, self._downsample)

    def _getSignal(self, mds_tree, tlim=None):
        'Get signal between `tlim` at sampling rate `self.Fs`.'
        # Retrieve full raw signal
        node = mds_tree.getNode(self._node_name)
        x = node.raw_of().data()

        # Determine time at the beginning of digitization record
        t0_dig = mds_tree.getNode(
            '.HARDWARE:%s:DI3' % self._digitizer_board).data()

        # Slice in time, if desired
        if (tlim is not None) or (self._downsample > 1):
            sl = self._getSlice(x, tlim=tlim, t0_dig=t0_dig)
            x = x[sl]
            t0 = t0_dig + (sl.start / (self.Fs * self._downsample))
        else:
            t0 = t0_dig

        return t0, x

    @property
    def t(self):
        pass
