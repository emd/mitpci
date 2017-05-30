from nose import tools
import numpy as np
from mitpci.interferometer.ellipse import FittedEllipse, _ellipse


def test__ellipse():
    E = np.arange(0, 2 * np.pi, np.pi / 180)

    # Standard unit circle
    x, y = _ellipse(1, 1, 0, 0, 0, E=E)
    np.testing.assert_array_equal(x, np.cos(E))
    np.testing.assert_array_equal(y, np.sin(E))

    # Offset unit circle
    x0 = 1.
    y0 = -2.5
    x, y = _ellipse(1, 1, x0, y0, 0, E=E)
    np.testing.assert_allclose(x, x0 + np.cos(E))
    np.testing.assert_allclose(y, y0 + np.sin(E))

    # Stretched
    a = 2
    b = 0.5
    x, y = _ellipse(a, b, 0, 0, 0, E=E)
    np.testing.assert_allclose(x, a * np.cos(E))
    np.testing.assert_allclose(y, b * np.sin(E))

    # Rotated unit circle
    x, y = _ellipse(1, 1, 0, 0, np.pi / 2, E=E)
    np.testing.assert_allclose(x, -np.sin(E), atol=1e-16)
    np.testing.assert_allclose(y, np.cos(E), atol=1e-16)

    return


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
