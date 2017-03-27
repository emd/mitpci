'''This module implements a class for correlating the MIT heterodyne
interferometer (toroidal location: 285 degrees, R = 1.98 m) with the
DIII-D V2 heterodyne interferometer (toroidal location: 240 degrees,
R = 1.94 m).

'''


# Standard library imports
# import numpy as np
# import matplotlib.pyplot as plt

# Related 3rd-party imports
from .demodulated import Demodulated
import bci
from random_data.spectra import CrossSpectralDensity


class ToroidalCorrelation(CrossSpectralDensity):
    def __init__(self, D, V2=None, trigger_offset=None,
                 vibration_subtracted=True):
        self._checkInputs(D, V2)

        # Load V2 data, if not provided by user
        if V2 is None:
            # Determine initial and final times
            # of MIT interferometer record
            tlim = D.I.t()[[0, -1]]

            # V2 record is *bracketed* by `tlim`; that is,
            # V2.t()[0] >= tlim[0] and V2.t()[-1] <= tlim[-1].
            V2 = bci.signal.Signal(
                D.shot, chord='V2', beam='CO2', tlim=tlim,
                vibration_subtracted=vibration_subtracted)

        self.shot = D.shot
        self.trigger_offset = trigger_offset
        self.vibration_subtracted = V2.vibration_subtracted

        # Interpolate MIT interferometer onto V2 time base

        # Compute cross-spectral density

        # Colormaps

    def _checkInputs(self, D, V2):
        'Check that `D` and `V2` are the correct types and are compatible.'
        # Valid types for `D` and `V2`
        Dtype = Demodulated
        V2type = bci.signal.Signal

        # Check that user-provided `D` is of correct type
        if type(D) is not Dtype:
            raise ValueError('`D` must be of type %s' % Dtype)

        # Check user-provided `V2` is of correct type and
        # is compatible with `D`
        if V2 is None:
            pass
        elif type(V2) is not V2type:
            raise ValueError('`V2` must be of type %s' % V2type)
        elif D.shot != V2.shot:
            raise ValueError('`D` and `V2` correspond to different shots')
        elif (D.I.t()[0] > V2.t()[-1]) or (D.I.t()[-1] < V2.t()[0]):
            raise ValueError('No temporal overlap between `D` and `V2`')

        return
