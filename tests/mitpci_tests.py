from nose import tools
import numpy as np
from mitpci.signal import Signal


def test_getDigitizerBoard():
    # Use the default "model" tree
    shot = -1

    # Boundaries of board 7
    tools.assert_equal(Signal(shot, 1)._digitizer_board, 'DT216_7')
    tools.assert_equal(Signal(shot, 8)._digitizer_board, 'DT216_7')

    # Boundaries of board 8
    tools.assert_equal(Signal(shot, 9)._digitizer_board, 'DT216_8')
    tools.assert_equal(Signal(shot, 16)._digitizer_board, 'DT216_8')

    # Other channels should raise ValueError
    tools.assert_raises(ValueError, Signal, *[shot, 0])
    tools.assert_raises(ValueError, Signal, *[shot, 17])


def test_getSampleRate():
    # Chris likes this shot - will probably be available for tests forever
    shot = 150000
    Fs = 4e6  # digitization rate

    # Without a specified sample rate, there will be no downsampling and
    # the sample rate will be equal to the intrinsic digitization rate
    sig = Signal(shot, 1)
    tools.assert_equal(sig.Fs, Fs)
    tools.assert_equal(sig._downsample, 1)

    # Integer downsampling specified
    downsample = 10
    sig = Signal(shot, 1, Fs=(Fs / downsample))
    tools.assert_equal(sig.Fs, Fs / downsample)
    tools.assert_equal(sig._downsample, downsample)

    # Non-integer downsampling specified; time resolution must be
    # at least as good as specified, so realized downsampling factor
    # will be the *floor* of the downsampling required to achieve
    # the specified `Fs`
    downsample = 2.5
    sig = Signal(shot, 1, Fs=(Fs / downsample))
    tools.assert_equal(sig.Fs, Fs / np.floor(downsample))
    tools.assert_equal(sig._downsample, np.floor(downsample))

    # "Negative" downsampling should raise a ValueError
    tools.assert_raises(ValueError, Signal, *[shot, 1], **{'Fs': -4e6})
