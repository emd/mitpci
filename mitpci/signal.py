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

        self._checkChannels(channels_per_board=channels_per_board)

        self._digitizer_board = self._getDigitizerBoard(
            channels_per_board=channels_per_board)
        # self._node_name = self.getNodeName(
        #     channels_per_board=channels_per_board)

        mds_tree = mds.Tree('pci', shot=shot, mode='ReadOnly')
        self.Fs, self._downsample = self._getSampleRate(mds_tree, Fs=Fs)

        pass

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

    # def _getNodeName(self, channels_per_board=8):
    #     board_channel = 1 + ((self.channel - 1) % channels_per_board)
    #     node_name = '.HARDWARE:%s:INPUT_%s' % (self._digitizer_board,
    #                                            str(board_channel).zfill(2))
    #     return node_name

    def getInitialTime(self):
        pass

    def getSignal(self):
        pass

    @property
    def t(self):
        pass
