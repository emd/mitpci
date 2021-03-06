from nose import tools
import numpy as np
import MDSplus as mds
from mitpci.signal import Signal


def test__checkChannel():
    # Use the default "model" tree
    shot = -1

    #  Channels < 0 or > 16 should raise ValueError
    tools.assert_raises(ValueError, Signal, *[shot, 0])
    tools.assert_raises(ValueError, Signal, *[shot, 17])

    return


def test__getDigitizerBoard():
    # Use the default "model" tree
    shot = -1

    # Boundaries of board 7
    tools.assert_equal(Signal(shot, 1)._digitizer_board, 'DT216_7')
    tools.assert_equal(Signal(shot, 8)._digitizer_board, 'DT216_7')

    # Boundaries of board 8
    tools.assert_equal(Signal(shot, 9)._digitizer_board, 'DT216_8')
    tools.assert_equal(Signal(shot, 16)._digitizer_board, 'DT216_8')

    return


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

    return


def test__getSampleRate():
    shot = 167340
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

    return


def test__getSlice():
    # Use the default "model" tree
    shot = -1

    # Test signal and limiting slicing values
    x = np.arange(10)
    imin = 0
    imax = len(x)

    # Create `Signal` object
    sig = Signal(shot, 1)

    # ValueError tests
    tools.assert_raises(ValueError, sig._getSlice, *[x], **{'tlim': [1]})
    tools.assert_raises(ValueError, sig._getSlice, *[x], **{'tlim': [1, 2, 3]})

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

    # -------------------------------------------------------------------------
    # (2) Downsampling, trivial `t0_dig`
    # -------------------------------------------------------------------------
    # Vary `sig.Fs` and `sig._downsample` inversely such that
    # the `tlim` values can stay the same as above.
    sig._downsample = 2  # Note: must be an *integer*
    sig.Fs /= sig._downsample

    # (2a) Downsampling only
    tools.assert_equal(
        sig._getSlice(x, tlim=None, t0_dig=0.),
        slice(imin, imax, sig._downsample))

    # (2b) Slicing from both ends, no downsampling
    tlim = [1, 8]
    tools.assert_equal(
        sig._getSlice(x, tlim=tlim, t0_dig=0.),
        slice(tlim[0], tlim[1] + 1, sig._downsample))

    # (2c) Slicing from lower end only, no downsampling
    tlim = [1, 10 + np.finfo(float).eps]
    tools.assert_equal(
        sig._getSlice(x, tlim=tlim, t0_dig=0.),
        slice(tlim[0], imax, sig._downsample))

    # (2d) Slicing from upper end only, no downsampling
    tlim = [-1 - np.finfo(float).eps, 8]
    tools.assert_equal(
        sig._getSlice(x, tlim=tlim, t0_dig=0.),
        slice(imin, tlim[1] + 1, sig._downsample))

    # -------------------------------------------------------------------------
    # (3) Downsampling, non-trivial `t0_dig`
    # -------------------------------------------------------------------------
    # Lower slicing indices *within* the bounds of the array should be
    # decremented by `np.floor(t0_dig)` relative to (2), and
    # upper slicing indices *within* the bounds of the array should be
    # decremented by `np.ceil(t0_dig)` (Note that this is equivalent to
    # (`1 + np.floor(t0_dig)`) relative to (2).
    t0_dig = 1.5

    # (3a) Downsampling only; `t0_dig` does not come into play
    # when `tlim` is None
    tools.assert_equal(
        sig._getSlice(x, tlim=None, t0_dig=t0_dig),
        slice(imin, imax, sig._downsample))

    # (3b) Slicing from both ends, no downsampling
    tlim = [2, 8]
    tools.assert_equal(
        sig._getSlice(x, tlim=tlim, t0_dig=t0_dig),
        slice(tlim[0] - np.floor(t0_dig),
              tlim[1] - np.floor(t0_dig),
              sig._downsample))

    # (3c) Slicing from lower end only, no downsampling
    tlim = [2, 10 + np.finfo(float).eps + t0_dig]
    tools.assert_equal(
        sig._getSlice(x, tlim=tlim, t0_dig=t0_dig),
        slice(tlim[0] - np.floor(t0_dig), imax, sig._downsample))

    # (3d) Slicing from upper end only, no downsampling
    tlim = [-1 - np.finfo(float).eps + t0_dig, 8]
    tools.assert_equal(
        sig._getSlice(x, tlim=tlim, t0_dig=t0_dig),
        slice(imin, tlim[1] - np.floor(t0_dig), sig._downsample))

    return


def test__getSignal():
    shot = 167340

    # Slice signal by (a) downsampling (Fs_dig = 4e6 in this shot) and
    # (b) looking at only a portion of the digitized record
    tlim = [1, 2]
    sig = Signal(shot, 1, Fs=1e6, tlim=tlim)

    # -------------------------------------------------------------------------
    # (1) Check `sig.t0`
    # -------------------------------------------------------------------------
    # `sig.t0` should be bracketed (inclusive) below by `tlim[0]` and
    # bracketed (non-inclusive) above by
    # `tlim[0] + (1. / (sig.Fs * sig._downsample))`
    tools.assert_less_equal(tlim[0], sig.t0)
    tools.assert_less(sig.t0, tlim[0] + (1. / (sig.Fs * sig._downsample)))

    # -------------------------------------------------------------------------
    # (2) Check that length of `sig.x` is correct
    # -------------------------------------------------------------------------
    Fs_dig = sig.Fs * sig._downsample  # digitization rate

    # `t_1` is the time corresponding to `x[-1]`
    t_1 = sig.t0 + np.floor((tlim[1] - sig.t0) * Fs_dig) / Fs_dig

    # Number of samples in *full* digitized record
    N_full = ((t_1 - sig.t0) * Fs_dig) + 1

    # Number of records in the retrieved, *downsampled* record
    N_retrieved = np.floor(N_full / sig._downsample)

    # For internal self-consistency. If `N_retrieved is *not* equal to
    # its integer representation, the above calculations for the expected
    # signal length are incorrect
    tools.assert_equal(
        N_retrieved,
        np.int(N_retrieved),
        msg='`N_retrieved` should be an integer; test calculations incorrect.')

    tools.assert_equal(N_retrieved, len(sig.x))

    return


def test_t():
    # Use the default "model" tree (Don't load signal from MDSplus)
    shot = -1

    # Create `Signal` object
    sig = Signal(shot, 1)

    # Now, overwrite properties of `sig` to easily test `_getSlice()` method
    sig.Fs = 1.
    sig._downsample = 1
    sig.t0 = 0
    sig.x = np.arange(10)

    # (1) `sig.t()` is equal to `sig.x` for the above parameters
    np.testing.assert_equal(sig.t(), sig.x)

    # (2) Linear shift
    sig.t0 += 1
    np.testing.assert_equal(sig.t(), sig.x + 1)

    # (3) Double the sampling rate at which signal is retrieved
    sig.Fs *= 2
    np.testing.assert_equal(sig.t(), (sig.x / sig.Fs) + 1)

    # (4) The time base computed by the `sig.t()` method corresponds
    # to the retrieved points in `sig.x`, which may be downsampled
    # from the full digitized record. The downsampling is folded
    # into the sampling rate `sig.Fs` at which the signal is retrieved,
    # however, so altering `sig._downsample` should *not* affect `sig.t()`;
    # that is, `sig.t()` should remain unchanged from (3)
    sig._downsample = 2
    np.testing.assert_equal(sig.t(), (sig.x / sig.Fs) + 1)

    return


def test_volts_per_bit():
    shot = 167340

    # Create `Signal` object
    sig = Signal(shot, 1)

    # Read bits-to-voltage conversion factor from tree.
    # Note that the conversion factor retrieved from the tree
    # applies uniformly to *all* of the channels on a given board,
    # while the `Signal.volts_per_bit` property allows specification
    # of the conversion factor on a channel-by-channel basis,
    # regardless of their board location.
    mds_tree = mds.Tree('pci', shot=shot, mode='ReadOnly')
    node = mds_tree.getNode('.HARDWARE:%s:TO_VOLTS' % sig._digitizer_board)
    tree_conversion_factor = node.data()

    tools.assert_equal(sig.volts_per_bit, tree_conversion_factor)

    return
