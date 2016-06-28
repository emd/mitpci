'''This module implements a class for retrieving, processing, and analyzing
signals from the interferometer channel of the mitpci system.

'''


# Standard library imports
import numpy as np
import matplotlib.pyplot as plt

# Related 3rd-party imports
from .signal import Signal


class Demodulated(object):
    '''An object corresponding to the I&Q signals of the interferometer.

    Attributes:
    -----------
    I - :py:class:`Signal <mitpci.signal.Signal>` instance
        The in-phase (I) signal of the interferometer.

    Q - :py:class:`Signal <mitpci.signal.Signal>` instance
        The quadrature (Q) signal of the interferometer.

    Methods:
    --------

    '''
    def __init__(
            self, shot, channel_I=1, channel_Q=2,
            quiet=False, **signal_kwargs):
        '''Create an instance of the `Demodulated` class.

        Input parameters:
        -----------------
        shot - int
            The shot number of the signal to be retrieved.

        channel_I [channel Q] - int
            The channel of the mitpci system corresponding to the
            interferometer's in-phase (I) [quadrature (Q) signal]signal

        quiet - bool
            If True, suppress printing messages to the terminal.

        signal_kwargs - any valid keyword arguments for
            :py:class:`Signal <mitpci.signal.Signal>`.

            For example, use

                    D = Demodulated(167340, tlim=[1, 3.3])

            to retrieve the I and Q signals for shot 167340 between
            1 <= t [s] <= 3.3 from channels 1 and 2, respectively,
            of the mitpci system.

        '''
        self.shot = shot

        if not quiet:
            print '\nRetrieving in-phase (I) signal for %i' % shot
        self.I = Signal(shot, channel_I, **signal_kwargs)

        if not quiet:
            print 'Retrieving quadrature (Q) signal for %i' % shot
        self.Q = Signal(shot, channel_Q, **signal_kwargs)

        return
