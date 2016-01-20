from nose import tools
from mitpci.signal import Signal


def test_getDigitizerBoard():
    # Use the default "model" tree
    shot = -1

    # Boundaries of board 7
    tools.assert_equal(Signal(shot, 1).digitizer_board, 'DT216_7')
    tools.assert_equal(Signal(shot, 8).digitizer_board, 'DT216_7')

    # Boundaries of board 8
    tools.assert_equal(Signal(shot, 9).digitizer_board, 'DT216_8')
    tools.assert_equal(Signal(shot, 16).digitizer_board, 'DT216_8')

    # Other channels should raise ValueError
    tools.assert_raises(ValueError, Signal(shot, 0))
    tools.assert_raises(ValueError, Signal(shot, 16))
