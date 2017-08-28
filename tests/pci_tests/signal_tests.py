from nose import tools
import numpy as np
import filters
from mitpci.pci.signal import Phase


def test_Phase__init__():
    shot = 171110
    channel = [8]

    # Needs to receive `filt` that is `None` or `filters.fir.Kaiser`
    tools.assert_raises(
        ValueError,
        Phase, *[shot, channel], **{'filt': {}})

    # Use incorrect sampling rate
    Fs_bad = 1e6
    hpf = filters.fir.Kaiser(-120, 5e3, 10e3, pass_zero=False, Fs=Fs_bad)
    tools.assert_raises(
        ValueError,
        Phase, *[shot, channel], **{'filt': hpf})

    # Check that multi-channel retrieval is equivalent to multiple
    # single-channel retrievals and that the records are arranged
    # as expected
    tlim = [1, 2]
    Ph7 = Phase(shot, 7, tlim=tlim)
    Ph8 = Phase(shot, 8, tlim=tlim)
    Ph = Phase(shot, [7, 8], tlim=tlim)

    np.testing.assert_equal(Ph7.x, Ph.x[0, :])
    np.testing.assert_equal(Ph8.x, Ph.x[1, :])

    return
