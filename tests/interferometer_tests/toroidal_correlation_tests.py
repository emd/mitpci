import numpy as np
from nose import tools
from mitpci.interferometer import Demodulated, ToroidalCorrelation
import bci


def test__checkSignals():
    shot = 167341
    tlim = np.array([2., 2.25])
    dt = tlim[1] - tlim[0]

    # Load data
    D = Demodulated(shot, tlim=tlim)
    V2 = bci.signal.Signal(shot, tlim=tlim)

    # Need to pass object of type :py:class:`Demodulated
    # <mitpci.signal.demodulated.Demodulated>`, not just signal
    tools.assert_raises(ValueError, ToroidalCorrelation, D.getPhase())

    # Need to pass :py:class:`Signal <bci.signal.Signal>` for V2 data
    tools.assert_raises(ValueError, ToroidalCorrelation, D, {'V2': V2.x})

    # `D` and `V2` must correspond to same shot
    V2.shot += 1
    tools.assert_raises(ValueError, ToroidalCorrelation, D, {'V2': V2})

    # `V2` temporal range just below that of `D`
    tlim_lo = tlim - dt - (1. / V2.Fs)
    V2 = bci.signal.Signal(shot, tlim=tlim_lo)
    tools.assert_raises(ValueError, ToroidalCorrelation, D, {'V2': V2})

    # `V2` temporal range just above that of `D`
    tlim_hi = tlim + dt + (1. / V2.Fs)
    V2 = bci.signal.Signal(shot, tlim=tlim_hi)
    tools.assert_raises(ValueError, ToroidalCorrelation, D, {'V2': V2})

    return
