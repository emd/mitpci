'''This module implements a class for determining the mapping from
the MITPCI digitizer channels to the PCI detector array.

'''


# Standard library imports
import numpy as np
import matplotlib.pyplot as plt
import MDSplus as mds

# Related 3rd-party imports
from random_data.array import ArrayStencil


class DetectorArray(ArrayStencil):
    def __init__(self, shot, digitizer_channels):
        '''Create instance of the PCI `DetectorArray` class.

        Input parameters:
        -----------------
        shot - int
            The DIII-D shot number.

        digitizer_channels - array_like, (`N`,)
            The desired channels of the MITPCI digitizer system to include
            in the analysis of PCI measurements. This may be the full set
            of available channels or any subset thereof. Digitizer channels
            corresponding to dead detector elements or mapped to something
            other than the PCI detector array (e.g. the I&Q signals of the
            heterodyne interferometer) are automatically removed. A
            ValueError is raised if non-existent digitizer channels are
            specified.

        '''
        self.shot = shot

        # Get digitizer to detector mapping and parse results
        mapping = self.getDigitizerToDetectorMapping(digitizer_channels)
        self.digitizer_channels = mapping[0]
        self.detector_elements = mapping[1]

        # # Create corresponding stencil
        # ArrayStencil.__init__(
        #     self,
        #     digitizer_channels,
        #     include_autocorrelations=True)

    def getDigitizerToDetectorMapping(self, digitizer_channels):
        # Determine the detector elements in use for `self.shot`.
        # Note that `detector_elements[i]` gives the detector
        # element corresponding to digitizer channel `i + 1`.
        tree = mds.Tree('pci', shot=self.shot, mode='ReadOnly')
        node = tree.getNode('ELECTRONICS:DETECTOR_CH')
        detector_elements = node.getData().data()

        min_digitizer_channel = 1
        max_digitizer_channel = len(detector_elements)

        # Ensure that values in `digitizer_channels` are valid & sorted
        if ((np.min(digitizer_channels) < min_digitizer_channel) or
                np.max(digitizer_channels) > max_digitizer_channel):
            raise ValueError(
                'Values in `digitizer_channels` must be between %i and %i' %
                (min_digitizer_channel, max_digitizer_channel))
        else:
            digitizer_channels = np.sort(digitizer_channels)

        # Only consider detector elements corresponding to
        # the requested digitizer channels, remembering that
        # `detector_elements[i]` corresponds to the `i + 1`
        # digitizer channel.
        detector_elements = detector_elements[digitizer_channels - 1]

        # Remove digitizer channels corresponding to dead detector elements
        # or mapped to something other than the PCI detector array (e.g.
        # the I&Q signals of the heterodyne interferometer)
        ind = np.where(detector_elements != 0)[0]
        detector_elements = detector_elements[ind]
        digitizer_channels = digitizer_channels[ind]

        return digitizer_channels, detector_elements
