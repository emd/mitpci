'''This module implements a class for fitting a collection
of points (x, y) to an ellipse. If x = x(t) and y = y(t) and
the ellipse evolves adiabatically (i.e. much slower than
x and y traverse the ellipse), then the ellipse's dynamic
parameters can be extracted.

'''


# # Standard library imports
import numpy as np
# import matplotlib.pyplot as plt

# Related 3rd-party imports
from fit_ellipse import fit_ellipse


class FittedEllipse(object):
    def __init__(self, x, y, starts):
        '''Create an instance of the `FittedEllipse` class.

        Parameters:
        -----------
        x - array_like, (`N`,)
            The x coordinate of the points to be fit to an ellipse.
            [x] = AU

        y - array_like, (`N`,)
            The y coordinate of the points to be fit to an ellipse.
            [y] = AU

        starts - array_like, (`M`,)
            The points described by `x` and `y` are subdivided
            into `M` smaller segments, with the i-th segment
            obtained via slicing with

                    slice(starts[i], starts[i + 1])

            The points in the i-th segment are then fit
            to determine the properties of the i-th ellipse.
            In total, `M` ellipses will be fit.

        '''
        # Note the *unique* identity of the `x` and `y` arrays
        self.xid = id(x)
        self.yid = id(y)

        # Ellipse fitting
        a, b, x0, y0, theta = self.getFits(x, y, starts)
        self.a = a
        self.b = b
        self.x0 = x0
        self.y0 = y0
        self.theta = theta

    def getSlice(self, ellipse, starts):
        'Get slice for ellipse w/ initial index `starts[ellipse]`.'
        if ellipse != (len(starts) - 1):
            sl = slice(starts[ellipse], starts[ellipse + 1])
        else:
            # The final ellipse is a boundary case that
            # requires special handling
            sl = slice(starts[ellipse], None)

        return sl

    def getFits(self, x, y, starts):
        'Get fits for ellipses w/ initial indices `starts`.'
        # Initialize arrays
        N = len(starts)
        a = np.zeros(N)
        b = np.zeros(N)
        x0 = np.zeros(N)
        y0 = np.zeros(N)
        theta = np.zeros(N)

        # Fit each ellipse
        for e in np.arange(N):
            sl = self.getSlice(e, starts)
            a[e], b[e], x0[e], y0[e], theta[e]= fit_ellipse(x[sl], y[sl])

        return a, b, x0, y0, theta


def _ellipse(a, b, x0, y0, theta, E=np.arange(0, 2 * np.pi, np.pi / 180)):
    '''Get (x, y) coordinates of specified ellipse.

    Parameters:
    -----------
    a - float
        Semi-major axis of ellipse.
        [a] = AU

    b - float
        Semi-minor axis of ellipse.
        [b] = [a]

    x0 - float
        x-coordinate of ellipse center.
        [x0] = [a]

    y0 - float
        y-coordinate of ellipse center.
        [y0] = [a]

    theta - float
        Angle between semi-major axis of ellipse and x-axis.
        [theta] = rad

    E - array_like, (`N`,)
        The "eccentric anomaly" of the ellipse.
        [E] = rad

    Returns:
    --------
    (x, y) - tuple, where:

    x - array_like, (`N`,)
        x-coordinates of ellipse

    y - array_like, (`N`,)
        y-coordinates of ellipse

    '''
    # Build up x-coordinate iteratively
    x = a * np.cos(theta) * np.cos(E)
    x -= (b * np.sin(theta) * np.sin(E))
    x += x0

    # Build up y-coordinate iteratively
    y = a * np.sin(theta) * np.cos(E)
    y += (b * np.cos(theta) * np.sin(E))
    y += y0

    return x, y
