'''This module implements a class for fitting a collection
of points (x, y) to an ellipse. If x = x(t) and y = y(t) and
the ellipse evolves adiabatically (i.e. much slower than
x and y traverse the ellipse), then the ellipse's dynamic
parameters can be extracted.

'''


# # Standard library imports
import numpy as np
import matplotlib.pyplot as plt

# Related 3rd-party imports
from fit_ellipse import fit_ellipse


class FittedEllipse(object):
    def __init__(self, x, y, starts, t=None):
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

        t - array_like, (`N`,)
            The timestamps corresponding to points in `x` and `y`.
            The timestamp separation is assumed to be constant.
            [t] = AU

        '''
        if len(x) != len(y):
            raise ValueError('`x` and `y` must have same length')

        # Note the *unique* identity of the `x` and `y` arrays
        self.xid = id(x)
        self.yid = id(y)

        self.starts = starts

        # Ellipse fitting
        a, b, x0, y0, theta = self.getFits(x, y)
        self.a = a
        self.b = b
        self.x0 = x0
        self.y0 = y0
        self.theta = theta

        # Determine times corresponding to midpoint of each fitted slice
        if t is not None:
            if len(t) == len(x):
                self.t = self.getSliceTimes(t)
            else:
                raise ValueError('`t` must have same length as `x` and `y`')
        else:
            self.t = None

    def getSlice(self, ellipse):
        'Get slice for ellipse w/ initial index `self.starts[ellipse]`.'
        if ellipse != (len(self.starts) - 1):
            sl = slice(self.starts[ellipse], self.starts[ellipse + 1])
        else:
            # The final ellipse is a boundary case that
            # requires special handling
            sl = slice(self.starts[ellipse], None)

        return sl

    def getSliceTimes(self, t):
        'Get times corresponding to midpoint of each slice.'
        N = len(self.starts)
        tmid = np.zeros(N)

        for e in np.arange(N):
            sl = self.getSlice(e)
            tmid[e] = 0.5 * (t[sl][-1] + t[sl][0])

        return tmid

    def getFits(self, x, y):
        'Get fits for ellipses w/ initial indices `self.starts`.'
        # Initialize arrays
        N = len(self.starts)
        a = np.zeros(N)
        b = np.zeros(N)
        x0 = np.zeros(N)
        y0 = np.zeros(N)
        theta = np.zeros(N)

        # Fit each ellipse
        for e in np.arange(N):
            sl = self.getSlice(e)
            a[e], b[e], x0[e], y0[e], theta[e] = fit_ellipse(x[sl], y[sl])

        return a, b, x0, y0, theta

    def plotFits(self, tlim=None):
        fig, axes = plt.subplots(5, 1, sharex=True)

        # Create a local variable for dependent axis
        if self.t is not None:
            t = self.t.copy()
        else:
            t = np.arange(len(self.a))

        if tlim is not None:
            tind = np.arange(np.logical_and(
                t >= tlim[0],
                t <= tlim[-1]))[0]
        else:
            tind = slice(None, None)

        axes[0].plot(t[tind], self.a[tind])
        axes[0].set_ylabel(r'$a$')

        axes[1].plot(t[tind], self.b[tind])
        axes[1].set_ylabel(r'$b$')

        axes[2].plot(t[tind], self.x0[tind])
        axes[2].set_ylabel(r'$x_0$')

        axes[3].plot(t[tind], self.y0[tind])
        axes[3].set_ylabel(r'$y_0$')

        axes[4].plot(t[tind], self.theta[tind])
        axes[4].set_ylabel(r'$\theta$')
        axes[4].set_xlabel(r'$t$')

        plt.show()

        return fig, axes

    def compensateEllipticity(self, x, y):
        '''Use fitting parameters to invert affine transformation,
        effectively mapping the ellipse back to a circle.

        '''
        if id(x) != self.xid:
            raise ValueError('Provided `x` *not* equal to fitted `x`')
        if id(y) != self.yid:
            raise ValueError('Provided `y` *not* equal to fitted `y`')

        # Initialize compensated arrays
        xc = np.zeros(len(x))
        yc = np.zeros(len(y))

        # Map each ellipse back to its corresponding circle
        for e in np.arange(len(self.starts)):
            sl = self.getSlice(e)

            # Radius of corresponding circle is simply
            # an average of semi-major and semi-minor axes
            r = 0.5 * (self.a[e] + self.b[e])

            # Build up compensated x-coordinate iteratively
            xc[sl] = (x[sl] - self.x0[e]) * np.cos(self.theta[e])
            xc[sl] += ((y[sl] - self.y0[e]) * np.sin(self.theta[e]))
            xc[sl] *= (r / self.a[e])

            # Build up compensated y-coordinate iteratively
            yc[sl] = (y[sl] - self.y0[e]) * np.cos(self.theta[e])
            yc[sl] -= ((x[sl] - self.x0[e]) * np.sin(self.theta[e]))
            yc[sl] *= (r / self.b[e])

        return xc, yc


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
