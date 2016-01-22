from nose import tools
import numpy as np
from mitpci.signal import Signal


def test__checkChannels():
    # Use the default "model" tree
    shot = -1

    #  Channels < 0 or > 16 should raise ValueError
    tools.assert_raises(ValueError, Signal, *[shot, 0])
    tools.assert_raises(ValueError, Signal, *[shot, 17])


def test__getDigitizerBoard():
    # Use the default "model" tree
    shot = -1

    # Boundaries of board 7
    tools.assert_equal(Signal(shot, 1)._digitizer_board, 'DT216_7')
    tools.assert_equal(Signal(shot, 8)._digitizer_board, 'DT216_7')

    # Boundaries of board 8
    tools.assert_equal(Signal(shot, 9)._digitizer_board, 'DT216_8')
    tools.assert_equal(Signal(shot, 16)._digitizer_board, 'DT216_8')


def test__getNodeName():
    # Use the default "model" tree
    shot = -1

    # Boundaries of board 7
    tools.assert_equal(
        Signal(shot, 1)._node_name,
        '.HARDWARE:DT216_7:INPUT_01')
    tools.assert_equal(
        Signal(shot, 8)._node_name,
        '.HARDWARE:DT216_7:INPUT_08')

    # Boundaries of board 8
    tools.assert_equal(
        Signal(shot, 9)._node_name,
        '.HARDWARE:DT216_8:INPUT_01')
    tools.assert_equal(
        Signal(shot, 16)._node_name,
        '.HARDWARE:DT216_8:INPUT_08')


def test__getSampleRate():
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


def test_getSlice():
    # Use the default "model" tree
    shot = -1

    # Test signal and limiting slicing values
    x = np.arange(10)
    imin = 0
    imax = len(x)

    # Create `Signal` object
    sig = Signal(shot, 1)

    # -------------------------------------------------------------------------
    # (1) No downsampling, trivial `t0_dig`
    # -------------------------------------------------------------------------
    # Now, overwrite properties of `sig` to easily test `_getSlice()` method
    sig.Fs = 1.
    sig._downsample = 1

    # (1a) No actual slicing
    tools.assert_equal(
        sig._getSlice(x, tlim=None, t0_dig=0.),
        slice(imin, imax, 1))

    # (1b) Slicing from both ends, no downsampling
    tlim = [1, 8]
    tools.assert_equal(
        sig._getSlice(x, tlim=tlim, t0_dig=0.),
        slice(tlim[0], tlim[1] + 1, 1))

    # (1c) Slicing from lower end only, no downsampling
    tlim = [1, 10 + np.finfo(float).eps]
    tools.assert_equal(
        sig._getSlice(x, tlim=tlim, t0_dig=0.),
        slice(tlim[0], imax, 1))

    # (1d) Slicing from upper end only, no downsampling
    tlim = [-1 - np.finfo(float).eps, 8]
    tools.assert_equal(
        sig._getSlice(x, tlim=tlim, t0_dig=0.),
        slice(imin, tlim[1] + 1, 1))

    # # -------------------------------------------------------------------------
    # # (2) Downsampling, trivial `t0_dig`
    # # -------------------------------------------------------------------------
    # # Vary `sig.Fs` and `sig._downsample` inversely such that
    # # the `tlim` values can stay the same as above.
    # sig._downsample = 2  # Note: must be an *integer*
    # sig.Fs /= sig._downsample

    # # (2a) Downsampling only
    # tools.assert_equal(
    #     sig._getSlice(x, tlim=None, t0_dig=0.),
    #     slice(0, len(x) + 1, sig._downsample))

    # # (2b) Slicing from both ends, downsampling
    # tlim = [1, 9]
    # tools.assert_equal(
    #     sig._getSlice(x, tlim=tlim, t0_dig=0.),
    #     slice(tlim[0], tlim[1] + 1, sig._downsample))

    # # (2c) Slicing from lower end only, downsampling
    # tlim = [1, 10.1]
    # tools.assert_equal(
    #     sig._getSlice(x, tlim=tlim, t0_dig=0.),
    #     slice(tlim[0], imax, sig._downsample))

    # # (2d) Slicing from upper end only, downsampling
    # tlim = [-0.1, 9]
    # tools.assert_equal(
    #     sig._getSlice(x, tlim=tlim, t0_dig=0.),
    #     slice(imin, tlim[1], sig._downsample))

    # # -------------------------------------------------------------------------
    # # (3) Downsampling, non-trivial `t0_dig`
    # # -------------------------------------------------------------------------
    # # Note: `np.ceil(t0_dig) = 1` implies that all slice start and stop values
    # # should be shiftedj
    # t0_dig = 0.5

    # # (3a) Downsampling only; `t0_dig` should *not* affect this computation,
    # # as `t0_dig` only comes into play when `tlim` is not None
    # tools.assert_equal(
    #     sig._getSlice(x, tlim=None, t0_dig=t0_dig),
    #     slice(0, len(x) + 1, sig._downsample))

    # # (3b) Slicing from both ends, downsampling
    # tlim = [1, 9]
    # tools.assert_equal(
    #     sig._getSlice(x, tlim=tlim, t0_dig=t0_dig),
    #     slice(tlim[0], tlim[1] + 1, sig._downsample))

    # # (3c) Slicing from lower end only, downsampling
    # tlim = [1, 10.1]
    # tools.assert_equal(
    #     sig._getSlice(x, tlim=tlim, t0_dig=t0_dig),
    #     slice(tlim[0], imax, sig._downsample))

    # # (3d) Slicing from upper end only, downsampling
    # tlim = [-0.1, 9]
    # tools.assert_equal(
    #     sig._getSlice(x, tlim=tlim, t0_dig=t0_dig),
    #     slice(imin, tlim[1], sig._downsample))


def test__getSignal():
    # downsampling - how to test easily/effectively?
    # value error for wrong size t_lim
    # slicing limits - for both inside and outside record length,
    #   on both sides of record; how to test easily/effectively?
    # returned t0 value
    pass
