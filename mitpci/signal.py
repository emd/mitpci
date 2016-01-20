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
    def __init__(self, shot, channel=None, Fs=None, tlim=None):
        pass
