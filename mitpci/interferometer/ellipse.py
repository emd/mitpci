'''This module implements a class for fitting a collection
of points (x, y) to an ellipse. If x = x(t) and y = y(t) and
the ellipse evolves adiabatically (i.e. much slower than
x and y traverse the ellipse), then the ellipse's dynamic
parameters can be extracted.

'''


# # Standard library imports
# import numpy as np
# import matplotlib.pyplot as plt

# Related 3rd-party imports
from fit_ellipse import fit_ellipse


class FittedEllipse(object):
    def __init__(self, x, y, starts):
        pass

    def getSlice(self, ellipse, starts):
        'Get slice for ellipse w/ initial index `starts[ellipse]`.'
        if ellipse != (len(starts) - 1):
            sl = slice(starts[ellipse], starts[ellipse + 1])
        else:
            # The final ellipse is a boundary case that
            # requires special handling
            sl = slice(starts[ellipse], None)

        return sl

    def getFit(self, x, y, starts):
        # for sl in slices:
        #     a, b, center0, center1, phi = fit_ellipse(x, y)
        #     center, axes = (center0, center1), (a, b)

        return
