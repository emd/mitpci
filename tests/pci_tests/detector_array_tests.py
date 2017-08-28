from nose import tools
import numpy as np
import MDSplus as mds
from mitpci.pci.detector_array import DetectorArray


def test_DetectorArray__init__():
    # Ensure value errors are raised when invalid digitizer channels
    # are specified
    shot = -1
    tree = mds.Tree('pci', shot=shot, mode='ReadOnly')
    node = tree.getNode('ELECTRONICS:DETECTOR_CH')
    detector_elements = node.getData().data()
    max_digitizer_channel = len(detector_elements)

    # Digitizer channel too low
    tools.assert_raises(
        ValueError,
        DetectorArray,
        *(shot, [0]))

    # Digitizer channel too high
    tools.assert_raises(
        ValueError,
        DetectorArray,
        *(shot, [max_digitizer_channel + 1]))

    # Ensure that we get correct mapping for a past shot
    shot = 171521

    # Detector elements, raw mapping straight from MDSplus:
    #
    #   [0, 0, 2, 7, 10, 13, 15, 16, 17, 0, 0, 18, 20, 23, 26, 31]
    #
    # But we know that zeros should not be included in retrieved mapping.
    # The expected mapping from digitizer channels to PCI detector
    # elements is:
    digitizer_channels_exp = [3, 4,  5,  6,  7,  8,  9, 12, 13, 14, 15, 16]
    detector_elements_exp =  [2, 7, 10, 13, 15, 16, 17, 18, 20, 23, 26, 31]

    # Let's use a "naive" list of all digitizer channels (including
    # dead channels and channels not mapped to the PCI detector), and
    # let's additionally make the list non-monotonic. If `DetectorArray`
    # is working correctly, we should still get agreement with
    # the expected mapping listed above.
    #
    # First, create non-monotonic digitizer-channel listing.
    digitizer_channels_naive = np.arange(max_digitizer_channel) + 1
    ind = len(digitizer_channels_naive) // 2
    digitizer_channels_naive = np.concatenate((
        digitizer_channels_naive[ind:],
        digitizer_channels_naive[0:ind]))

    # Build detector from naive digitizer-channel listing
    detector = DetectorArray(shot, digitizer_channels_naive)

    np.testing.assert_equal(
        detector.digitizer_channels,
        digitizer_channels_exp)

    np.testing.assert_equal(
        detector.detector_elements,
        detector_elements_exp)

    return
