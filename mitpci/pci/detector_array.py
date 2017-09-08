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
    '''An object characterizing the PCI detector-to-digitizer mapping.

    This class is derived from :py:class:`ArrayStencil
    <random_data.array.ArrayStencil>` and thus shares all of its
    attributes and methods. Attributes and properties *unique* to
    :py:class:`DetectorArray <mitpci.pci.detector_array.DetectorArray>`
    are discussed below.

    Attributes:
    -----------
    shot - int
        The DIII-D shot number.

    digitizer_channels - array_like, (`M`,)
        The channels of the MITPCI digitizer system that map to the
        detector elements in `self.locations`.

    Properties:
    -----------
    detector_elements - array_like, (`M`,)
        A property that simply points to `self.locations`, which gives
        the detector elements that correspond to `self.digitizer_channels`.

    '''
    def __init__(self, shot, digitizer_channels,
                 include_autocorrelations=True):
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

        include_autocorrelations - bool
            If True, include autocorrelations as unique correlation pairs.

        '''
        self.shot = shot

        # Get digitizer to detector mapping and parse results. Note that
        # `digitizer_channels` is explicitly a class attribute, while
        # `detector_elements` is not. This is because the subsequent
        # call to the base class's (i.e. `ArrayStencil`) initialization
        # routine creates an attribute called `locations` that is
        # identical to `detector_elements`; to make this interdependence
        # explicit, we've defined a class property `detector_elements()`
        # that simply points to `locations`.
        mapping = self.getDigitizerToDetectorMapping(digitizer_channels)
        self.digitizer_channels = mapping[0]
        detector_elements = mapping[1]

        # Create corresponding stencil
        ArrayStencil.__init__(
            self,
            detector_elements,
            include_autocorrelations=include_autocorrelations)

    def getDigitizerToDetectorMapping(self, digitizer_channels):
        'Get mapping of digitizer channels to PCI detector elements.'
        # Determine the detector elements in use for `self.shot`.
        # Note that `detector_elements[i]` gives the detector
        # element corresponding to digitizer channel `i + 1`.
        tree = mds.Tree('pci', shot=self.shot, mode='ReadOnly')
        node = tree.getNode('ELECTRONICS:DETECTOR_CH')
        detector_elements = node.getData().data()

        min_digitizer_channel = 1
        max_digitizer_channel = len(detector_elements)

        # Ensure that values in `digitizer_channels` are valid & sorted.
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
        # the I&Q signals of the heterodyne interferometer).
        ind = np.where(detector_elements != 0)[0]
        detector_elements = detector_elements[ind]
        digitizer_channels = digitizer_channels[ind]

        # Finally, sort `detector_elements` from smallest to largest;
        # apply the same indexing to `digitizer_channels` to preserve
        # the correct detector-to-digitizer mapping. By sorting here,
        # we prevent any subsequent chronological sorting of
        # `detector_elements` from re-ordering `detector_elements`
        # relative to `digitizer_channels`, which would destroy the
        # detector-to-digitizer mapping.
        ind = np.argsort(detector_elements)
        detector_elements = detector_elements[ind]
        digitizer_channels = digitizer_channels[ind]

        return digitizer_channels, detector_elements

    @property
    def detector_elements(self):
        'Get detector elements corresponding to `self.digitizer_channels.'
        return self.locations
