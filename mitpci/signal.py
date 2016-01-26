'''This module implements a class for retrieving signals digitized by
the mitpci system.

'''


import numpy as np
import MDSplus as mds


class Signal(object):
    '''An object corresponding to the signal retrieved from the mitpci system.

    Attributes:
    -----------
    shot - int
        The shot number of the retrieved signal.

    channel - int
        The channel number of the retrieved signal.

    x - array-like, (`N`,)
        The retrieved signal. Note that this is the "raw" signal and
        has the units of bits.
        [x] = bits

    Fs - float
        The sampling rate at which the signal was retrieved. Note that this
        may be smaller than the signal's initial digitization rate; that is,
        the retrieved signal may be downsampled relative to the full signal.
        [Fs] = samples / second

    t0 - float
        The time corresponding to the first retrieved point in the signal;
        that is, if x(t) corresponds to the continuous signal being sampled,
        then `Signal.x[0]` = x(t0)
        [t0] = s

    t - array-like, (`N`,)
        The times corresponding to the points in the retrieved signal.
        This attribute is actually implemented as a class property
        so that the time-base is generated on the fly as needed and
        is not stored as an object property; this helps save memory and
        processing time.
        [t] = s

    '''
    def __init__(self, shot, channel, channels_per_board=8,
                 Fs=None, tlim=None):
        '''Create an instance of the `Signal` class.

        Input parameters:
        -----------------
        shot - int
            The shot number of the signal to be retrieved.

        channel - int
            The channel number of the signal to be retrieved.
            Currently, there are two digitizer boards, each of which
            digitizes `channels_per_board` channels. A ValueError is
            raised if `channel <= 0` or `channel > 2 * channels_per_board`.

        channels_per_board - int
           The number of channels digitized per digitizer board.
           Currently, there are two digitizer boards, each of which
           digitizes 8 channels. The value of `channels_per_board`
           determines whether the value of `channel` is valid or not.

        Fs - float
            The desired sampling rate at which to retrieve the signal.
            This may differ from the original digitization rate, `Fs_dig`;
            that is, the retrieved signal can be downsampled. While
            the specified `Fs` can be any positive value, the internal
            methods of this class will process the user's request and
            produce a *realized* sampling rate `Signal.Fs` that is
            subject to the following two constraints:

                (1) `Signal.Fs <= Fs_dig`, and
                (2) `Signal. Fs / Fs_dig = m`, where `m` is an integer

            `Fs <= 0` raises a ValueError.

            Note: the "raw" signal is returned from the MDSplus server.
            Josh Stillerman and Tom Fredian (of MIT, and the developers
            of MDSplus) informed me that there is no way to read in only
            a portion of the raw data. Thus, the entire raw signal must be
            read into memory, and the signal will subsequently be sliced
            in time if needed. As a result, specifying a small `Fs` will
            not save computational time but will end up saving memory.

            [Fs] = samples / second

        tlim - array-like, (2,)
            The lower and upper limits in time for which the signal
            will be retrieved.

            The specified `tlim` values will always bracket the retrieved
            signal. That is, if `tlim[0]` does not correspond to an exact
            digitization time, then the initial time returned (`Signal.t0`)
            will be the closest digitization time *greater* than `tlim[0]`.
            Similarly, if `tlim[1]` does not correspond to an exact
            digitization time, then the final time (`Signal.t[-1]`) will be
            the closest digitization time *less* than `tlim[1]`. Further,
            if `tlim[0]` is less than the initial digitization time,
            the retrieved signal will begin with the initial digitized point.
            Similarly, if `tlim[1]` exceeds the final digitization time,
            the retrieved signal will end with the final digitized point.

            A ValueError is raised if `len(tlim) != 2`.

            Note: the "raw" signal is returned from the MDSplus server.
            Josh Stillerman and Tom Fredian (of MIT, and the developers
            of MDSplus) informed me that there is no way to read in only
            a portion of the raw data. Thus, the entire raw signal must be
            read into memory, and the signal will subsequently be sliced
            in time if needed. As a result, specifying a small `tlim` will
            not save computational time but will end up saving memory.

            [tlim] = s

        '''
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

    def t(self):
        'Get times for points in `self.x`.'
        return self.t0 + (np.arange(len(self.x)) / self.Fs)
