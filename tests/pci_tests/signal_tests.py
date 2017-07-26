from nose import tools
import filters
from mitpci.pci.signal import Phase


def test_Phase__init__():
    shot = 171110
    channel = 8

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

    return
