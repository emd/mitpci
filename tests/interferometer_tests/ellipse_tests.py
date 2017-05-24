from nose import tools
import numpy as np
from mitpci.interferometer.ellipse import FittedEllipse


def test_FittedEllipse_getSlice():
    # Unit circle w/ 1-degree spacing between successive points
    # and 2 full transits around origin
    N = 2
    dth = 1
    th = np.arange(0, 360 * N, dth)
    x = np.cos((np.pi / 180) * th)
    y = np.sin((np.pi / 180) * th)
    starts = np.where((th % 360) == 0)[0]

    E = FittedEllipse(x, y, starts)

    # Check slicing for first ellipse
    tools.assert_equal(
        E.getSlice(0, starts),
        slice(0, 360))

    # Check slicing for final ellipse
    tools.assert_equal(
        E.getSlice(1, starts),
        slice(360, None))

    return
