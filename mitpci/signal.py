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
    def __init__(self, shot, channel, Fs=None, tlim=None):
        self.shot = shot
        self.channel = channel

        self.digitizer_board = self.getDigitizerBoard()

        # mds_tree = mds.Tree('pci', shot=shot, mode='ReadOnly')
        pass

    def getDigitizerBoard(self, channels_per_board=8):
        'Get digitizer board corresponding to `self.channel`.'
        if (self.channel > 0) and (self.channel <= channels_per_board):
            return 'DT216_7'
        elif (self.channel > 0) and (self.channel <= (2 * channels_per_board)):
            return 'DT216_8'
        else:
            raise ValueError(
                'Valid `channel` values between 1 and %i' %
                (2 * channels_per_board))

    def getSampleRate(self, mds_tree):
        # need digitizer...
        # Fs = mds_tree.getNode(self._digitizers[0] + ':CLOCK_DIV').data()
        pass

    def getInitialTime(self):
        pass

    def getSignal(self):
        pass
