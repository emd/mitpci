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
    # and 1 full transits around origin
    N = 1
    dE = 1
    E = np.arange(0, 360 * N, dE)
    x, y = _ellipse(1, 1, 0, 0, 0, E=E)
    starts = np.where((E % 360) == 0)[0]

    FE = FittedEllipse(x, y, starts)

    # Check slicing for first (and only) ellipse
    tools.assert_equal(
        FE.getSlice(0, starts),
        slice(0, None))

    # Unit circle w/ 1-degree spacing between successive points
    # and 2 full transits around origin
    N = 2
    dE = 1
    E = np.arange(0, 360 * N, dE)
    x, y = _ellipse(1, 1, 0, 0, 0, E=E)
    starts = np.where((E % 360) == 0)[0]

    FE = FittedEllipse(x, y, starts)

    # Check slicing for first ellipse
    tools.assert_equal(
        FE.getSlice(0, starts),
        slice(0, 360))

    # Check slicing for final ellipse
    tools.assert_equal(
        FE.getSlice(1, starts),
        slice(360, None))

    return


def test_FittedEllipse_getFits():
    E = np.arange(0, 2 * np.pi, np.pi / 180)

    # Single ellipse
    # Note: Not quite sure on angle convention
    # of ellipse-fitting algorithm, so just use
    # small, positive value (which ellipse-fitting
    # algorithm maps back onto the same value)
    axes1 = np.abs(np.random.randn(2))
    a1 = np.max(axes1)
    b1 = np.min(axes1)
    x01 = np.random.randn()
    y01 = np.random.randn()
    theta1 = 0.1
    x1, y1 = _ellipse(a1, b1, x01, y01, theta1, E=E)

    FE = FittedEllipse(x1, y1, [0])

    np.testing.assert_allclose(FE.a, np.array([a1]))
    np.testing.assert_allclose(FE.b, np.array([b1]))
    np.testing.assert_allclose(FE.x0, np.array([x01]))
    np.testing.assert_allclose(FE.y0, np.array([y01]))
    np.testing.assert_allclose(
        FE.theta,
        np.array([theta1]))

    # Generate 2nd, distinct ellipse
    # Note: Not quite sure on angle convention
    # of ellipse-fitting algorithm, so just use
    # small, positive value (which ellipse-fitting
    # algorithm maps back onto the same value)
    axes2 = np.abs(np.random.randn(2))
    a2 = np.max(axes2)
    b2 = np.min(axes2)
    x02 = np.random.randn()
    y02 = np.random.randn()
    theta2 = 0.2
    x2, y2 = _ellipse(a2, b2, x02, y02, theta2, E=E)

    FE = FittedEllipse(
        np.concatenate((x1, x2)),
        np.concatenate((y1, y2)),
        [0, len(E)])

    np.testing.assert_allclose(FE.a, np.array([a1, a2]))
    np.testing.assert_allclose(FE.b, np.array([b1, b2]))
    np.testing.assert_allclose(FE.x0, np.array([x01, x02]))
    np.testing.assert_allclose(FE.y0, np.array([y01, y02]))
    np.testing.assert_allclose(FE.theta, np.array([theta1, theta2]))

    return


def test_FittedEllipse_compensateEllipticity_WrongInputs():
    E = np.arange(0, 2 * np.pi, np.pi / 180)

    # Single ellipse
    # Note: Not quite sure on angle convention
    # of ellipse-fitting algorithm, so just use
    # small, positive value (which ellipse-fitting
    # algorithm maps back onto the same value)
    axes1 = np.abs(np.random.randn(2))
    a1 = np.max(axes1)
    b1 = np.min(axes1)
    x01 = np.random.randn()
    y01 = np.random.randn()
    theta1 = 0.1
    x1, y1 = _ellipse(a1, b1, x01, y01, theta1, E=E)
    starts = [0]

    FE = FittedEllipse(x1, y1, starts)

    tools.assert_raises(
        ValueError,
        FE.compensateEllipticity, *[x1, x1, starts])

    tools.assert_raises(
        ValueError,
        FE.compensateEllipticity, *[y1, y1, starts])

    return


def test_FittedEllipse_compensateEllipticity():
    # Irrational angular spacing to avoid zero crossings
    # in sine and cosine functions, which throw off
    # the relative-tolerance algorithms below
    # (i.e. eps / 0 -> infty, which is not small)
    dE = np.sqrt(2) * (np.pi / 180)
    E = np.arange(dE, 2 * np.pi, dE)

    # Single ellipse
    # Note: Not quite sure on angle convention
    # of ellipse-fitting algorithm, so just use
    # small, positive value (which ellipse-fitting
    # algorithm maps back onto the same value)
    axes1 = np.abs(np.random.randn(2))
    a1 = np.max(axes1)
    b1 = np.min(axes1)
    x01 = np.random.randn()
    y01 = np.random.randn()
    theta1 = 0.1
    x1, y1 = _ellipse(a1, b1, x01, y01, theta1, E=E)

    # Fit elliptical data and apply compensation
    starts = [0]
    FE = FittedEllipse(x1, y1, starts)
    xc, yc = FE.compensateEllipticity(x1, y1, starts)

    # Compare to expectations
    r1 = 0.5 * (a1 + b1)
    np.testing.assert_allclose(xc, r1 * np.cos(E))
    np.testing.assert_allclose(yc, r1 * np.sin(E))

    # Generate 2nd, distinct ellipse
    # Note: Not quite sure on angle convention
    # of ellipse-fitting algorithm, so just use
    # small, positive value (which ellipse-fitting
    # algorithm maps back onto the same value)
    axes2 = np.abs(np.random.randn(2))
    a2 = np.max(axes2)
    b2 = np.min(axes2)
    x02 = np.random.randn()
    y02 = np.random.randn()
    theta2 = 0.2
    x2, y2 = _ellipse(a2, b2, x02, y02, theta2, E=E)

    # Fit elliptical data and apply compensation
    x = np.concatenate((x1, x2))  # need unique identifier/memory
    y = np.concatenate((y1, y2))  # need unique identifier/memory
    starts = [0, len(E)]
    FE = FittedEllipse(x, y, starts)
    xc, yc = FE.compensateEllipticity(x, y, starts)

    # Compare to expectations
    sl1 = FE.getSlice(0, starts)
    np.testing.assert_allclose(xc[sl1], r1 * np.cos(E))
    np.testing.assert_allclose(yc[sl1], r1 * np.sin(E))

    r2 = 0.5 * (a2 + b2)
    sl2 = FE.getSlice(1, starts)
    np.testing.assert_allclose(xc[sl2], r2 * np.cos(E))
    np.testing.assert_allclose(yc[sl2], r2 * np.sin(E))

    return
