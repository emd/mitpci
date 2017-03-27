'''This module implements a class for correlating the MIT heterodyne
interferometer (toroidal location: 285 degrees, R = 1.98 m) with the
DIII-D V2 heterodyne interferometer (toroidal location: 240 degrees,
R = 1.94 m).

'''


# Standard library imports
import numpy as np
from scipy.interpolate import interp1d

# Related 3rd-party imports
from .demodulated import Demodulated
import bci
from random_data.spectra import CrossSpectralDensity


class ToroidalCorrelation(CrossSpectralDensity):
    '''A class for toroidally correlating the MIT and V2 interferometers.

    This class is derived from :py:class:`CrossSpectralDensity
    <random_data.spectra.CrossSpectralDensity>` and thus shares
    all of its attributes and methods. Attributes *unique*
    to :py:class:`ToroidalCorrelation
    <mitpci.interferometer.toroidal_correlation.ToroidalCorrelation>`
    are discussed below.

    Attributes:
    -----------
    shot - int
        DIII-D shot number.

    trigger_offset - float
        The timebase offset between the MIT and V2 interferometers.
        (See `__init__()` docstring for more information).

    vibration_subtracted - bool
        If True, vibrational contributions from the V2-measured phase
        have been removed.

    '''
    def __init__(self, D, V2=None, trigger_offset=-32.5e-6,
                 vibration_subtracted=True, **csd_kwargs):
        '''Create an instance of the `ToroidalCorrelation` class.

        Input Parameters:
        -----------------
        D - :py:class:`Demodulated
                <mitpci.interferometer.demodulated.Demodulated>`
            Object corresponding to the MIT interferometer's
            demodulated in-phase (I) and quadrature (Q) signals.

        V2 - None, or :py:class:`Signal <bci.signal.Signal>`
            Object corresponding to V2 interferometer's phase.
            If `V2` is `None`, the V2-measured phase corresponding
            to `D` is automatically loaded.

        trigger_offset - float
            The timebase offset between the MIT and V2 interferometers.
            Although the MIT and V2 interferometer clocks are phase-locked,
            discrepancies between nominal and realized trigger times and
            nominal and actual sampling rates in both systems can result
            in a finite "trigger offset". This offset must be compensated
            for in order to extract sensible results from the correlation.

            When `trigger_offset` is positive, the MIT timebase leads
            the V2 timebase. When `trigger_offset` is negative, the
            MIT timebase lags the V2 timebase.

            [trigger_offset] = s

        vibration_subtracted - bool
            If True, remove vibrational contributions from the
            V2-measured phase. Vibrational contributions to the
            CO2-measured phase are typically "small" for frequencies
            above 10 kHz. While removing the vibrational contributions
            from the CO2-measured phase can increase the signal-to-noise
            ratio in *some cases*, it can potentially introduce other
            problems, as is discussed on the BCI homepage.

        csd_kwargs - any valid keyword arguments for
            :py:class:`CrossSpectralDensity
                <random_data.spectra.CrossSpectralDensity>`.

            For example, use

                    xcorr = ToroidalCorrelation(...,
                        Tens=5e-3, Nreal_per_ens=10)

            to specify 5-ms ensembles with 10 realizations per ensemble.
            See the `CrossSpectralDensity` documentation for
            further details.

            Note that the `t0` and `Fs` keywords will be
            neglected if provided in `csd_kwargs`, as these
            parameters are automatically determined from
            `D` and `V2`.

        '''
        self._checkSignals(D, V2)
        csd_kwargs = self._checkCsdKwargs(csd_kwargs)

        # Load V2 data, if not provided by user
        if V2 is None:
            # Determine initial and final times
            # of MIT interferometer record
            tlim = D.I.t()[[0, -1]]

            # V2 record is *bracketed* by `tlim`; that is,
            # V2.t()[0] >= tlim[0] and V2.t()[-1] <= tlim[-1].
            V2 = bci.signal.Signal(
                D.shot, chord='V2', beam='CO2', tlim=tlim,
                vibration_subtracted=vibration_subtracted)

        self.shot = D.shot
        self.trigger_offset = trigger_offset
        self.vibration_subtracted = V2.vibration_subtracted

        # Interpolate MIT interferometer onto V2 time base
        print '\nInterpolating MIT measurements onto V2 timebase'
        ph_MIT_interp = (interp1d(D.I.t(), D.getPhase()))(V2.t())

        # Account for trigger offset
        sl_V2, sl_MIT = self._getOffsetSlices(V2)

        # Compute cross-spectral density
        CrossSpectralDensity.__init__(
            self, V2.x[sl_V2], ph_MIT_interp[sl_MIT],
            Fs=V2.Fs, t0=V2.t0,
            **csd_kwargs)

        # Colormaps and easier plotting

    def _checkSignals(self, D, V2):
        'Check that `D` and `V2` are the correct types and are compatible.'
        # Valid types for `D` and `V2`
        Dtype = Demodulated
        V2type = bci.signal.Signal

        # Check that user-provided `D` is of correct type
        if type(D) is not Dtype:
            raise ValueError('`D` must be of type %s' % Dtype)

        # Check user-provided `V2` is of correct type and
        # is compatible with `D`
        if V2 is None:
            pass
        elif type(V2) is not V2type:
            raise ValueError('`V2` must be of type %s' % V2type)
        elif D.shot != V2.shot:
            raise ValueError('`D` and `V2` correspond to different shots')
        elif (D.I.t()[0] > V2.t()[-1]) or (D.I.t()[-1] < V2.t()[0]):
            raise ValueError('No temporal overlap between `D` and `V2`')

        return

    def _checkCsdKwargs(self, csd_kwargs):
        'Remove `t0` and `Fs` from `csd_kwargs`, if present.'
        forbidden_keys = ['t0', 'Fs']

        for key in forbidden_keys:
            if csd_kwargs.has_key(key):
                print '\nUser-specified `%s` is being neglected' % key
                csd_kwargs.pop(key)

        return csd_kwargs

    def _getOffsetSlices(self, V2):
        '''Get slices to shift V2 and MIT signals relative to one another
        in order to eliminate finite trigger offset between the two systems.

        '''
        offset = np.int(V2.Fs * self.trigger_offset)

        if offset > 0:
            sl_V2 = slice(None, -offset)
            sl_MIT = slice(offset, None)
        elif offset < 0:
            sl_V2 = slice(-offset, None)
            sl_MIT = slice(None, offset)
        else:
            sl_V2 = slice(None, None)
            sl_MIT = slice(None, None)

        return sl_V2, sl_MIT

    def plotModeNumber(self, nmin=-3, **plot_kwargs):
        '''Plot toroidal mode number as a function of frequency and time.

        Parameters:
        -----------
        nmin - int
            The minimum toroidal mode number to be plotted.
            The 45-degree toroidal separation of the V2 and MIT
            interferometers allows identification of 8 distinct
            toroidal mode numbers. The exact mode-number range,
            however, depends on the mode's rotation:

            - for unknown rotation, choose `nmin` = -3, which
              will yield a maximum toroidal mode number of 4
              (alternatively, `nmin` = -4 is equally sensible,
              and this yields a maximum toroidal mode number of 3);

            - for positive rotation (counterclockwise when viewing
              the vacuum vessel from above), choose `nmin` = 0,
              which yields a maximum toroidal mode number of 7;

        plot_kwargs - any valid keyword arguments for
            :py:instancemethod:`plotPhaseAngle
            <random_data.spectra.CrossSpectralDensity.plotPhaseAngle>`.

            Note that the {`dtheta`, `theta_min`, `mode_number`}
            keywords are prescribed in this method, so user-provided
            values for these keywords are overwritten.

        '''
        # Toroidal spacing of V2 and MIT interferometers
        # [dzeta] = radians
        dzeta = np.pi / 4

        # Phase angle corresponding to `nmin`
        theta_min = nmin * dzeta

        plot_kwargs['dtheta'] = dzeta
        plot_kwargs['theta_min'] = theta_min
        plot_kwargs['mode_number'] = True

        ax = CrossSpectralDensity.plotPhaseAngle(self, **plot_kwargs)

        return ax
